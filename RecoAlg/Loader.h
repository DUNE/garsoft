#include <string>
#include <vector>

#include "TFile.h"

std::string pnfs2xrootd(std::string loc, bool unauth = false);
/// Find files matching a UNIX glob, plus expand environment variables
std::vector<std::string> Wildcard(const std::string& wildcardString);

/// \brief Interface class for accessing ROOT files in sequence
class IFileSource
{
    public:
    virtual ~IFileSource() {}
    /// \brief Returns the next file in sequence, ready for reading
    ///
    /// A null return means that the end of the sequence has been reached.
    /// DO NOT close or delete the file that is returned.
    virtual TFile* GetNextFile() = 0;

    /// May return -1 indicating the number of files is not known
    virtual int NFiles() const {return -1;}
};

/// Simple file source based on an explicit list provided by the user
class FileListSource: public IFileSource
{
    public:
        FileListSource(const std::vector<std::string>& files);
        ~FileListSource();

        TFile* GetNextFile() override;
        int NFiles() const override {return fFileNames.size();}

        const std::vector<std::string>& GetFileNames() const { return fFileNames; }

    protected:
        std::vector<std::string> fFileNames; ///< The list of files
        std::vector<std::string>::iterator fIt; ///< Iterator into \ref fFileNames
        std::vector<std::string> fRetry; ///< List of files that failed 1st attempt
        bool fInRetry; ///< Did we finish fFileNames and are now in fRetry?
        TFile* fFile; ///< The most-recently-returned file
        static bool fgGotTickets; ///< Have we renewed our tickets?
};

// File source based on a wildcard (glob)
class WildcardSource: public FileListSource
{
    public:
        /// Wildcard or glob. Anything glob() accepts is OK
        /// May be a single literal filename
        WildcardSource(const std::string& wildcard);
        ~WildcardSource();

    protected:
        std::vector<std::string> CheckedWildcard(const std::string& wildcard) const;
};