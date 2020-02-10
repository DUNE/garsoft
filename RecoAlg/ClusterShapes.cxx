#include "ClusterShapes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

namespace util{

    ClusterShapes::ClusterShapes(int nhits, float* a, float* x, float* y, float* z)
    : _nHits(nhits), _aHit(nhits, 0.0), _xHit(nhits, 0.0), _yHit(nhits, 0.0), _zHit(nhits, 0.0), _ifNotGravity(1), _ifNotWidth(1), _ifNotInertia(1), _ifNotElipsoid(1)
    {
        for (int i(0); i < nhits; ++i)
        {
            _aHit[i] = a[i];
            _tHit[i] = 0.;
            _xHit[i] = x[i];
            _yHit[i] = y[i];
            _zHit[i] = z[i];
        }
    }

    ClusterShapes::ClusterShapes(int nhits, float* a, float* t, float* x, float* y, float* z)
    : _nHits(nhits), _aHit(nhits, 0.0), _tHit(nhits, 0.0), _xHit(nhits, 0.0), _yHit(nhits, 0.0), _zHit(nhits, 0.0), _ifNotGravity(1), _ifNotWidth(1), _ifNotInertia(1), _ifNotElipsoid(1)
    {
        for (int i(0); i < nhits; ++i)
        {
            _aHit[i] = a[i];
            _tHit[i] = t[i];
            _xHit[i] = x[i];
            _yHit[i] = y[i];
            _zHit[i] = z[i];
        }
    }

    //----------------------------------------------------------------------------

    ClusterShapes::ClusterShapes(int nhits, std::vector<float> a, std::vector<float> x, std::vector<float> y, std::vector<float> z)
    : _nHits(nhits), _aHit(nhits, 0.0), _xHit(nhits, 0.0), _yHit(nhits, 0.0), _zHit(nhits, 0.0), _ifNotGravity(1), _ifNotWidth(1), _ifNotInertia(1), _ifNotElipsoid(1)
    {
        for (int i(0); i < nhits; ++i)
        {
            _aHit[i] = a[i];
            _tHit[i] = 0.;
            _xHit[i] = x[i];
            _yHit[i] = y[i];
            _zHit[i] = z[i];
        }
    }

    ClusterShapes::ClusterShapes(int nhits, std::vector<float> a, std::vector<float> t, std::vector<float> x, std::vector<float> y, std::vector<float> z)
    : _nHits(nhits), _aHit(nhits, 0.0), _tHit(nhits, 0.0), _xHit(nhits, 0.0), _yHit(nhits, 0.0), _zHit(nhits, 0.0), _ifNotGravity(1), _ifNotWidth(1), _ifNotInertia(1), _ifNotElipsoid(1)
    {
        for (int i(0); i < nhits; ++i)
        {
            _aHit[i] = a[i];
            _tHit[i] = t[i];
            _xHit[i] = x[i];
            _yHit[i] = y[i];
            _zHit[i] = z[i];
        }
    }


    //----------------------------------------------------------------------------

    ClusterShapes::~ClusterShapes() {

    }

    //----------------------------------------------------------------------------

    int ClusterShapes::getNumberOfHits() {
        return _nHits;
    }

    //----------------------------------------------------------------------------

    float ClusterShapes::getTotalAmplitude() {
        if (_ifNotGravity == 1) findGravity();
        return _totAmpl;
    }

    //----------------------------------------------------------------------------

    float ClusterShapes::getAverageTime() {
        if (_ifNotGravity == 1) findGravity();
        return _totTime/_nHits;
    }

    //----------------------------------------------------------------------------

    float* ClusterShapes::getCenterOfGravity() {
        if (_ifNotGravity == 1) findGravity() ;
        return &_analogGravity[0] ;
    }

    //----------------------------------------------------------------------------

    float* ClusterShapes::getEigenValInertia() {
        if (_ifNotInertia == 1) findInertia();
        return &_ValAnalogInertia[0] ;
    }

    //----------------------------------------------------------------------------

    float* ClusterShapes::getEigenVecInertia() {
        if (_ifNotInertia == 1) findInertia();
        return &_VecAnalogInertia[0] ;
    }

    //----------------------------------------------------------------------------

    float ClusterShapes::getWidth() {
        if (_ifNotWidth == 1) findWidth();
        return _analogWidth;
    }

    void ClusterShapes::findElipsoid() {

        /**   Elipsoid parameter calculations see cluster_proper.f  */
        float cx,cy,cz ;
        float dx,dy,dz ;
        float r_hit_max, d_begn, d_last, r_max, proj;

        if (_ifNotInertia == 1) findInertia() ;

        //   Normalize the eigen values of inertia tensor
        float wr1 = sqrt(_ValAnalogInertia[0]/_totAmpl);
        float wr2 = sqrt(_ValAnalogInertia[1]/_totAmpl);
        float wr3 = sqrt(_ValAnalogInertia[2]/_totAmpl);

        _r1 = sqrt(wr2*wr3);                // spatial axis length -- the largest
        _r2 = sqrt(wr1*wr3);                // spatial axis length -- less
        _r3 = sqrt(wr1*wr2);                // spatial axis length -- even more less
        _vol = 4.*M_PI*_r1*_r2*_r3/3.;      // ellipsoid volume
        _r_ave = pow(_vol,1/3);             // average radius  (quibc root)
        _density = _totAmpl/_vol;           // density
        _eccentricity =_analogWidth/_r1;   // Cluster Eccentricity

        // Find Minumal and Maximal Lenght for Principal axis
        r_hit_max = -100000.;
        d_begn    =  100000.;
        d_last    = -100000.;
        cx = _VecAnalogInertia[0] ;
        cy = _VecAnalogInertia[1] ;
        cz = _VecAnalogInertia[2] ;

        for (int i(0); i < _nHits; ++i)
        {
            dx = _xHit[i] - _xgr;
            dy = _yHit[i] - _ygr;
            dz = _zHit[i] - _zgr;
            r_max = sqrt(dx*dx + dy*dy + dz*dz);

            if(r_max > r_hit_max) r_hit_max = r_max;
            proj = dx*cx + dy*cy + dz*cz;

            if(proj < d_begn)
            d_begn = proj;

            if(proj > d_last)
            d_last = proj;
        }

        _r1_forw = fabs(d_last);
        _r1_back = fabs(d_begn);

        _ifNotElipsoid = 0;
    }

    //----------------------------------------------------------------------------

    void ClusterShapes::findGravity() {

        _totAmpl = 0. ;

        for (int i(0); i < 3; ++i)
        _analogGravity[i] = 0.0 ;

        for (int i(0); i < _nHits; ++i)
        {
            _totAmpl+=_aHit[i] ;
            _totTime+=_tHit[i] ;
            _analogGravity[0]+=_aHit[i]*_xHit[i] ;
            _analogGravity[1]+=_aHit[i]*_yHit[i] ;
            _analogGravity[2]+=_aHit[i]*_zHit[i] ;
        }

        for (int i(0); i < 3; ++i)
        _analogGravity[i]/=_totAmpl ;

        _xgr = _analogGravity[0];
        _ygr = _analogGravity[1];
        _zgr = _analogGravity[2];

        _ifNotGravity = 0;
    }

    //----------------------------------------------------------------------------

    void ClusterShapes::findInertia() {

        double aIne[3][3];
        float radius2 = 0.0;

        findGravity();

        for (int i(0); i < 3; ++i) {
            for (int j(0); j < 3; ++j) {
                aIne[i][j] = 0.0;
            }
        }

        for (int i(0); i < _nHits; ++i)
        {
            float dX = _xHit[i] - _analogGravity[0];
            float dY = _yHit[i] - _analogGravity[1];
            float dZ = _zHit[i] - _analogGravity[2];
            aIne[0][0] += _aHit[i]*(dY*dY+dZ*dZ);
            aIne[1][1] += _aHit[i]*(dX*dX+dZ*dZ);
            aIne[2][2] += _aHit[i]*(dX*dX+dY*dY);
            aIne[0][1] -= _aHit[i]*dX*dY;
            aIne[0][2] -= _aHit[i]*dX*dZ;
            aIne[1][2] -= _aHit[i]*dY*dZ;
        }

        for (int i(0); i < 2; ++i) {
            for (int j = i+1; j < 3; ++j) {
                aIne[j][i] = aIne[i][j];
            }
        }
        //****************************************
        // analog Inertia
        //****************************************

        gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)aIne,3,3);
        gsl_vector* aVector = gsl_vector_alloc(3);
        gsl_matrix* aEigenVec = gsl_matrix_alloc(3,3);
        gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);
        gsl_eigen_symmv(&aMatrix.matrix,aVector,aEigenVec,wa);
        gsl_eigen_symmv_free(wa);
        gsl_eigen_symmv_sort(aVector,aEigenVec,GSL_EIGEN_SORT_ABS_ASC);

        for (int i(0); i < 3; i++) {
            _ValAnalogInertia[i] = gsl_vector_get(aVector,i);
            for (int j(0); j < 3; j++) {
                _VecAnalogInertia[i+3*j] = gsl_matrix_get(aEigenVec,i,j);
            }
        }

        // Main principal points away from IP
        _radius = 0.;
        radius2 = 0.;

        for (int i(0); i < 3; ++i) {
            _radius += _analogGravity[i]*_analogGravity[i];
            radius2 += (_analogGravity[i]+_VecAnalogInertia[i])*(_analogGravity[i]+_VecAnalogInertia[i]);
        }

        if ( radius2 < _radius) {
            for (int i(0); i < 3; ++i)
            _VecAnalogInertia[i] = - _VecAnalogInertia[i];
        }

        _radius = sqrt(_radius);
        _ifNotInertia = 0;

        // The final job
        findWidth();
        findElipsoid();

        gsl_vector_free(aVector);
        gsl_matrix_free(aEigenVec);

    }

    //----------------------------------------------------------------------------

    void ClusterShapes::findWidth() {

        float dist = 0.0;

        if (_ifNotInertia == 1)  findInertia() ;

        _analogWidth  = 0.0 ;

        for (int i(0); i < _nHits; ++i)
        {
            dist = findDistance(i) ;
            _analogWidth+=_aHit[i]*dist*dist ;
        }

        _analogWidth  = sqrt(_analogWidth / _totAmpl) ;

        _ifNotWidth = 0 ;
    }

    //----------------------------------------------------------------------------

    float ClusterShapes::findDistance(int i) {

        float cx = 0.0;
        float cy = 0.0;
        float cz = 0.0;
        float dx = 0.0;
        float dy = 0.0;
        float dz = 0.0;

        cx = _VecAnalogInertia[0] ;
        cy = _VecAnalogInertia[1] ;
        cz = _VecAnalogInertia[2] ;
        dx = _analogGravity[0] - _xHit[i] ;
        dy = _analogGravity[1] - _yHit[i] ;
        dz = _analogGravity[2] - _zHit[i] ;

        float tx = cy*dz - cz*dy ;
        float ty = cz*dx - cx*dz ;
        float tz = cx*dy - cy*dx ;
        float tt = sqrt(tx*tx+ty*ty+tz*tz) ;
        float ti = sqrt(cx*cx+cy*cy+cz*cz) ;
        float f = tt / ti ;

        return f ;
    }
}
