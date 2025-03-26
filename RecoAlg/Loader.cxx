#include "Loader.h"

#include "sys/stat.h"
#include "wordexp.h"
#include <unistd.h>

#include <cassert>
#include <iostream>

std::vector<std::string> sorted(const std::vector<std::string>& v)
{
  std::vector<std::string> ret;
  ret = v;
  std::sort(ret.begin(), ret.end());
  return ret;
}

//----------------------------------------------------------------------
std::string pnfs2xrootd(std::string loc, bool unauth)
{
  static bool first = true;
  static bool onsite = false;

  if (first && unauth) {
    first = false;
    char chostname[255];
    gethostname(chostname, 255);
    std::string hostname = chostname;

    if (hostname.find("fnal.gov") != std::string::npos) {
      onsite = true;
      std::cout << "Using unauthenticated xrootd access (port 1095) while on-site, hostname: "
                << hostname << std::endl;
    }
    else {
      onsite = false;
      std::cout << "Using authenticated xrootd access (port 1094) access while off-site, hostname: "
                << hostname << std::endl;
    }
  }

  if (loc.rfind("/pnfs/", 0) == 0) { // ie begins with
    if (onsite && unauth)
      loc = std::string("root://fndca1.fnal.gov:1095//pnfs/fnal.gov/usr/") + &loc.c_str()[6];
    else
      loc = std::string("root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr/") + &loc.c_str()[6];
  }
  return loc;
}

//----------------------------------------------------------------------
std::vector<std::string> Wildcard(const std::string& wildcardString)
{
  // Expand environment variables and wildcards like the shell would
  wordexp_t p;
  const int status = wordexp(wildcardString.c_str(), &p, WRDE_SHOWERR);

  if (status != 0) {
    std::cerr << "Wildcard string '" << wildcardString << "' returned error " << status
              << " from wordexp()." << std::endl;
    return {};
  }

  std::vector<std::string> fileList;

  for (unsigned int i = 0; i < p.we_wordc; ++i) {
    // Check the file exists before adding it
    struct stat sb;
    if (stat(p.we_wordv[i], &sb) == 0) fileList.push_back(p.we_wordv[i]);
  }

  wordfree(&p);

  return fileList;
}

bool FileListSource::fgGotTickets = false;

//----------------------------------------------------------------------
FileListSource::FileListSource(const std::vector<std::string>& files) : fInRetry(false), fFile(0)
{

  for (unsigned int i = 0; i < files.size(); ++i) {
    fFileNames.push_back(files[i]);
  }

  fIt = fFileNames.begin();

  for (const std::string& loc : fFileNames) {
    if (loc.rfind("/pnfs/", 0) == 0) { // ie begins with
      if (!fgGotTickets) {
        // No kerberos ticket means no point trying to voms-proxy-init. It
        // likely also means we're in a grid job, where that would be
        // counterproductive anyway.
        if (system("klist -5 -s || klist -s") != 0) fgGotTickets = true;
      }

      if (!fgGotTickets) {
        // This comes from NovaGridUtils or duneutil
        system("setup_fnal_security -b");

        fgGotTickets = true;
        break;
      }
    }
  }
}

//----------------------------------------------------------------------
FileListSource::~FileListSource()
{
  delete fFile;
}

//----------------------------------------------------------------------
TFile* FileListSource::GetNextFile()
{
  // Tidy up the last file we gave, which the caller no longer needs
  delete fFile;
  fFile = 0;

  // Did we run out of files?
  if (fInRetry && fIt == fRetry.end()) return 0;
  if (fIt == fFileNames.end()) {
    if (fRetry.empty()) return 0;
    fIt = fRetry.begin();
    fInRetry = true;
  }

  // If the file is on pnfs rewrite it to an xrootd address
  std::string loc = *fIt;
  loc = pnfs2xrootd(loc); // no-op for non /pnfs locations

  if (fInRetry) std::cout << "Retrying " << loc << " which was previously deferred..." << std::endl;

  fFile = TFile::Open(loc.c_str()); // This pattern allows xrootd

  if (!fFile) {
    if (!fInRetry) {
      std::cout << "Unable to open " << loc << std::endl;
      std::cout << "Will skip this file and try again later." << std::endl;
      fRetry.push_back(*fIt);
      ++fIt;
      return GetNextFile();
    }
    else {
      const int Nretry = 3;
      std::cout << "Failed to open. Will retry " << Nretry << " more times" << std::endl;
      for (int i = 0; i < Nretry; ++i) {
        std::cout << "Attempt " << i << std::endl;
        fFile = TFile::Open(loc.c_str());
        if (fFile) {
          std::cout << "Success!" << std::endl;
          break;
        }
      } // end for i
      if (!fFile) {
        std::cout << "Still unable to read " << loc << std::endl;
        std::cout << "Aborting" << std::endl;
        abort();
      }
    }
  }

  // The above should have guaranteed this
  assert(fFile);

  ++fIt; // Move on to the next file, for the subsequent call

  return fFile;
}

//----------------------------------------------------------------------
WildcardSource::WildcardSource(const std::string& wildcard)
  : FileListSource(sorted(CheckedWildcard(wildcard)))
{}

//----------------------------------------------------------------------
WildcardSource::~WildcardSource() {}

//----------------------------------------------------------------------
std::vector<std::string> WildcardSource::CheckedWildcard(const std::string& wildcard) const
{
  std::vector<std::string> ret = Wildcard(wildcard);

  struct stat ss;
  // If we found nothing, it may be because pnfs isn't mounted.
  if (ret.empty() && wildcard.find("/pnfs/") == 0 && stat("/pnfs/", &ss) != 0) {

    std::cout << "No files matching " << wildcard
              << " but that's probably because /pnfs is not mounted on the current grid node. If "
                 "you have to use your own files, either try running it interactivly or try "
                 "sam_add_dataset to create one."
              << std::endl;
  }

  return ret;
}