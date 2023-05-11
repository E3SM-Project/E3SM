#include "bound_exchange.h"

void bound_exchange(real4d &f, int dimz, int i_1, int i_2, int j_1, int j_2, int id) {
  YAKL_SCOPE( ncrms  , ::ncrms);

  real1d buffer("buffer", (nx+ny)*3*nz*ncrms);
  int i1  = i_1-1;
  int i2  = i_2-1;
  int j1  = j_1-1;
  int j2  = j_2-1;
  int i1p = i1+1;
  int i2p = i2+1;
  int j1p = j1+1;
  int j2p = j2+1;
  int offx, offy;

  if        (id==1) {
    offx = offx_u;
    offy = offy_u;
  } else if (id==2) {
    offx = offx_v;
    offy = offy_v; 
  } else if (id==3) {
    offx = offx_w;
    offy = offy_w; 
  } else if (id==4) {
    offx = offx_s;
    offy = offy_s; 
  } else if (id==5) {
    offx = offx_d;
    offy = offy_d; 
  } else {
    std::cout << "Id set in bound_exchange incorrectly:" << std::endl;
    exit(-1);
  }

  if (RUN3D) {
    //"North" -> "South":
    //Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,i+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,i+offx,icrm) = buffer(n);
    });

    //"North-East" -> "South-West":
    //Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int iStart = nx-i1-1;
      int iInd = nx-i1-1+i;
      int iEnd = nx-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,iInd+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int iStart = -i1-1;
      int iInd = -i1-1+i;
      int iEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });

    // "South-East" -> "North-West":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=1-1; j<i+j2; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int iStart = nx-i1-1;
      int iInd = nx-i1-1+i;
      int iEnd = nx-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,iInd+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int iStart = -i1-1;
      int iInd = -i1-1+i;
      int iEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });

    // "South" -> "North":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<1+j2; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,i+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,i+offx,icrm) = buffer(n);
    });

    // "South-West" -> "North-East":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<1+j2; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int iStart = 1-1;
      int iInd = 1-1+i;
      int iEnd = 1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,iInd+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int iStart = nxp1-1;
      int iInd = nxp1-1+i;
      int iEnd = nxp1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });
    
    // To "North-West" -> "South-East":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int iStart = 1-1;
      int iInd = 1-1+i;
      int iEnd = 1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(k,jInd+offy,iInd+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int iStart = nxp1-1;
      int iInd = nxp1-1+i;
      int iEnd = nxp1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });
     
  } //Endif for run3d

  //  "East" -> "West":
  // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = nx-i1-1;
    int iInd = nx-i1-1+i;
    int iEnd = nx-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    buffer(n) = f(k,j+offy,iInd+offx,icrm);
  });
  
  //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = -i1-1;
    int iInd = -i1-1+i;
    int iEnd = 0-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    f(k,j+offy,iInd+offx,icrm) = buffer(n);
  });

  // "West" -> "East":
  // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = 1-1;
    int iInd = 1-1+i;
    int iEnd = 1+i2-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    buffer(n) = f(k,j+offy,iInd+offx,icrm);
  });

  //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = nxp1-1;
    int iInd = nxp1-1+i;
    int iEnd = nxp1+i2-1;

    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    f(k,j+offy,iInd+offx,icrm) = buffer(n);
  });

}

void bound_exchange(real5d &f, int offL,int dimz, int i_1, int i_2, int j_1, int j_2, int id) {
  YAKL_SCOPE( ncrms  , ::ncrms);

  real1d buffer("buffer", (nx+ny)*3*nz*ncrms);
  int i1  = i_1-1;
  int i2  = i_2-1;
  int j1  = j_1-1;
  int j2  = j_2-1;
  int i1p = i1+1;
  int i2p = i2+1;
  int j1p = j1+1;
  int j2p = j2+1;
  int offx, offy;

  if        (id==1) {
    offx = offx_u;
    offy = offy_u;
  } else if (id==2) {
    offx = offx_v;
    offy = offy_v; 
  } else if (id==3) {
    offx = offx_w;
    offy = offy_w; 
  } else if (id==4) {
    offx = offx_s;
    offy = offy_s; 
  } else if (id==5) {
    offx = offx_d;
    offy = offy_d; 
  } else {
    std::cout << "Id set in bound_exchange incorrectly:" << std::endl;
    exit(-1);
  }

  if (RUN3D) {
    //"North" -> "South":
    //Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,i+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,i+offx,icrm) = buffer(n);
    });

    //"North-East" -> "South-West":
    //Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int iStart = nx-i1-1;
      int iInd = nx-i1-1+i;
      int iEnd = nx-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,iInd+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int iStart = -i1-1;
      int iInd = -i1-1+i;
      int iEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });

    // "South-East" -> "North-West":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=1-1; j<i+j2; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int iStart = nx-i1-1;
      int iInd = nx-i1-1+i;
      int iEnd = nx-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,iInd+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int iStart = -i1-1;
      int iInd = -i1-1+i;
      int iEnd = 0-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });

    // "South" -> "North":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<1+j2; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,i+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int n = _IDX(0,ncrms-1,icrm, 0,nx-1,i, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,i+offx,icrm) = buffer(n);
    });

    // "South-West" -> "North-East":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<1+j2; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = 1-1;
      int jInd = 1-1+j;
      int jEnd = 1+j2-1;
      int iStart = 1-1;
      int iInd = 1-1+i;
      int iEnd = 1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,iInd+offx,icrm);
    });

    // Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=nyp1-1; j<nyp1+j2; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j2p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = nyp1-1;
      int jInd = nyp1-1+j;
      int jEnd = nyp1+j2-1;
      int iStart = nxp1-1;
      int iInd = nxp1-1+i;
      int iEnd = nxp1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });
    
    // To "North-West" -> "South-East":
    // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=ny-j1-1; j<ny; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = ny-j1-1;
      int jInd = ny-j1-1+j;
      int jEnd = ny-1;
      int iStart = 1-1;
      int iInd = 1-1+i;
      int iEnd = 1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      buffer(n) = f(offL,k,jInd+offy,iInd+offx,icrm);
    });

    //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=-j1-1; j<0; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(dimz,j1p,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int jStart = -j1-1;
      int jInd = -j1-1+j;
      int jEnd = 0-1;
      int iStart = nxp1-1;
      int iInd = nxp1-1+i;
      int iEnd = nxp1+i2-1;
      int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, jStart,jEnd,jInd, 0,dimz-1,k);
      f(offL,k,jInd+offy,iInd+offx,icrm) = buffer(n);
    });
     
  } //Endif for run3d

  //  "East" -> "West":
  // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=nx-i1-1; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = nx-i1-1;
    int iInd = nx-i1-1+i;
    int iEnd = nx-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    buffer(n) = f(offL,k,j+offy,iInd+offx,icrm);
  });
  
  //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=-i1-1; i<0; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i1p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = -i1-1;
    int iInd = -i1-1+i;
    int iEnd = 0-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    f(offL,k,j+offy,iInd+offx,icrm) = buffer(n);
  });

  // "West" -> "East":
  // Pack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = 1-1;
    int iInd = 1-1+i;
    int iEnd = 1+i2-1;
    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    buffer(n) = f(offL,k,j+offy,iInd+offx,icrm);
  });

  //Unpack
    // for (int k=0; k<dimz; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=nxp1-1; i<nxp1+i2; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(dimz,ny,i2p,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int iStart = nxp1-1;
    int iInd = nxp1-1+i;
    int iEnd = nxp1+i2-1;

    int n = _IDX(0,ncrms-1,icrm, iStart,iEnd,iInd, 0,ny-1,j, 0,dimz-1,k);
    f(offL,k,j+offy,iInd+offx,icrm) = buffer(n);
  });

}
