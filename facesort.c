#include "common.h"

/* 
* Assign wall faces to BDBox
*/
void assignFaceToBDBox(double *vX, double *vY, double *vZ, int *con1, int *con2, int *con3,
						int noOfVert, int noOfF, double domDX, double domDY, double domDZ, 
                  int xD, int yD, int type){
   for(int i=0; i<noOfF; i++){
      double centroidX = (vX[con1[i]] + vX[con2[i]] + vX[con3[i]])/3.;
      double centroidY = (vY[con1[i]] + vY[con2[i]] + vY[con3[i]])/3.;
      double centroidZ = (vZ[con1[i]] + vZ[con2[i]] + vZ[con3[i]])/3.;
      int iIndex = floor((centroidX-xmin)/domDX);
      int jIndex = floor((centroidY-ymin)/domDY);
      int kIndex = floor((centroidZ-zmin)/domDZ);
      int cI = iIndex + jIndex*xD + kIndex*xD*yD;
      // writeLog3Num("logfile10.log","domD ",domDX/lengthFactor,domDY/lengthFactor,domDZ/lengthFactor);
      // writeLog3Num("logfile10.log","xD, yD ",xD,yD,0);
      // writeLog3Num("logfile10.log","iIndex ",iIndex,jIndex,kIndex);
   
      if(type == 1){
         bdBox[cI].face[bdBox[cI].totalFaces] = i;
         bdBox[cI].totalFaces += 1;
         maxFaces = fmax(maxFaces,bdBox[cI].totalFaces);
      }
      // if(type == 2){
      //    bdBox[cI].faceR[bdBox[cI].totalFacesR] = i;
      //    bdBox[cI].totalFacesR += 1;
      //    maxRotFaces = fmax(maxRotFaces,bdBox[cI].totalFacesR);
      // }
   }
   if(type == 1){
      writeLogNum("wallmesh.log","initial MAX FACES ",maxFaces);
   }
   else{
      writeLogNum("wallmesh.log","initial MAX FACES ROTOR ",maxRotFaces);
   }
   refineBdBox(xD, yD, zDiv, domDX, domDY, domDZ, vX, vY, vZ, con1, con2, con3, type);
}

/* Go through all DEM cells and filter faces within the cell.
    Filtering is done by checking the distance from cell DEM center 
    to projection on face
*/

void refineBdBox(int xDv, int yDv, int zDv, double dDx, double dDy, double dDz,
   double *vX, double *vY, double *vZ, int *con1, int *con2, int *con3, int type){
   //When searching near faces each cell checks neighbouring cells
   //extent of neighbouring cells is decided by offDiv index
   //int offDiv;
   //double edgeLength;
   tbdBox = allocateBdBoxArray(xDv*yDv*zDv);
   for(int i=0; i<xDiv*yDiv*zDiv; i++){
      tbdBox[i].noOfParticles = 0;
      tbdBox[i].totalFaces = 0;
      tbdBox[i].totalFaces = 0;
      tbdBox[i].face = allocateIntArray(MAX_FACES);
        //tbdBox[i].faceR = allocateIntArray(MAX_FACES);
        //bdBox[i].totalFaces = 0;
   }
   if(type == 1){
      //edgeLength = maxEdgeLength;
      
      //Make a temporary copy of from bdBox to tbdBox
      for(int k=0; k<zDv; k++){
         for(int j=0; j<yDv; j++){
            for(int i=0; i<xDv; i++){
               int cI = i + j*xDv + k*xDv*yDv;
               tbdBox[cI].totalFaces = bdBox[cI].totalFaces;
               for(int m=0; m<bdBox[cI].totalFaces; m++){
                  tbdBox[cI].face[m] = bdBox[cI].face[m];
                  //tbdBox[cI].faceType[m] = bdBox[cI].faceType[m];
               }
               bdBox[cI].totalFaces = 0;
            }
         }
      }
      assignToBdBox(tbdBox, vX, vY, vZ, con1, con2, con3, 
            xDv, yDv, zDv, dDx, dDy, dDz, maxEdgeLength, type);       
   }
   // else{
   //    //edgeLength = maxRotEdgeLength;
   //    //Make a temporary copy of from bdBox to tbdBox
   //    for(int k=0; k<zDv; k++){
   //       for(int j=0; j<yDv; j++){
   //          for(int i=0; i<xDv; i++){
   //             int cI = i + j*xDv + k*xDv*yDv;
   //             tbdBox[cI].totalFaces = bdBox[cI].totalFacesR;
   //             for(int m=0; m<bdBox[cI].totalFacesR; m++){
   //                tbdBox[cI].face[m] = bdBox[cI].faceR[m];
   //                //tbdBox[cI].faceType[m] = wBox[cI].faceType[m];
   //             }
   //             bdBox[cI].totalFacesR = 0;
   //          }
   //       }
   //    }
   //    assignToBdBox(tbdBox, vX, vY, vZ, con1, con2, con3,
   //          xDv, yDv, zDv, dDx, dDy, dDz, maxRotEdgeLength, type);
   // }
   free(tbdBox);
}


void assignToBdBox(struct BdBox *bBox, double *vX, double *vY, double *vZ, 
               int *con1, int *con2, int *con3, 
               int xDv, int yDv, int zDv, double dDx, double dDy, double dDz,
               double edgeLength, int type){
    int offDiv = ceil(edgeLength*lengthFactor/dDx);
    offDiv = 1;

    //writeLogNum("wallmesh.log","maxEdgeLength ",maxEdgeLength);
    //writeLogNum("wallmesh.log","domainDx ",domainDx);
   writeLogNum("wallmesh.log","div ",offDiv);
   for(int k=0; k<zDv; k++){
      for(int j=0; j<yDv; j++){
         for(int i=0; i<xDv; i++){
            if(i>=0 && j>=0 && k>=0){
               if(i<xDv && j<yDv && k<zDv){
                  int cI = i + j*xDv + k*xDv*yDv;

                  double *cellCent = allocateDoubleArray(DIM);
                  double *dVec = allocateDoubleArray(DIM);
                  cellCent[0] = i*dDx + 0.5*dDx + xmin;
                  cellCent[1] = j*dDy + 0.5*dDy + ymin;
                  cellCent[2] = k*dDz + 0.5*dDz + zmin;

                  //Go through neighbouring cells for face hunting
                  for(int kk=k-offDiv; kk<k+offDiv+1; kk++){
                     for(int jj=j-offDiv; jj<j+offDiv+1; jj++){
                        for(int ii=i-offDiv; ii<i+offDiv+1; ii++){
                           if(ii>=0 && jj>=0 && kk>=0){
                              if(ii<xDv && jj<yDv && kk<zDv){
                                 int cII = ii + jj*xDv + kk*xDv*yDv;
                                            ///if(tbdBox[cII].totalFaces>0){
                                                            //}
                                 for(int m=0; m<bBox[cII].totalFaces; m++){  
                                                //writeLogNum("wallmesh.log","tbdBox[cII].totalFaces ",wF[tbdBox[cII].face[m]].node1[0]/lengthFactor);
                                    int fId = bBox[cII].face[m];
                                    double *nd1 = allocateDoubleArray(DIM);
                                    double *nd2 = allocateDoubleArray(DIM);
                                    double *nd3 = allocateDoubleArray(DIM);
                                    nd1[0] = vX[con1[fId]];
                                    nd1[1] = vY[con1[fId]];
                                    nd1[2] = vZ[con1[fId]];
                                    nd2[0] = vX[con2[fId]];
                                    nd2[1] = vY[con2[fId]];
                                    nd2[2] = vZ[con2[fId]];
                                    nd3[0] = vX[con3[fId]];
                                    nd3[1] = vY[con3[fId]];
                                    nd3[2] = vZ[con3[fId]];

                                    double *tempDist = allocateDoubleArray(DIM);
                                    tempDist[0] = cellCent[0]-((nd1[0]+nd2[0]+nd3[0])/3.);
                                    tempDist[1] = cellCent[1]-((nd1[1]+nd2[1]+nd3[1])/3.);
                                    tempDist[2] = cellCent[2]-((nd1[2]+nd2[2]+nd3[2])/3.);

                                    getProjection(cellCent, nd1, nd2, nd3,dVec);
                                    double dist = 0.5*sqrt(dDx*dDx+
                                                      dDy*dDy+dDz*dDz)+
                                                      0.5*largestParDia*lengthFactor;

                                    if(vecMag(dVec) < dist && vecMag(tempDist) < dist*1.5){
                                      addFaceToBDBox(bBox[cII].face[m], cI, type);
                                    }
                                    free(tempDist);
                                    free(nd1);
                                    free(nd2);
                                    free(nd3);
                                 }//end of total faces
                              }
                           }
                        }
                     }
                  }//end of neighbouring i,j,k
                  free(dVec);
                  free(cellCent);
                  if(type == 1){
                     maxFaces = fmax(maxFaces,bdBox[cI].totalFaces); 
                  }
                  else{
                     maxRotFaces = fmax(maxRotFaces,bdBox[cI].totalFacesR); 
                  }
               }//end of i<iDiv && j<yDiv && k<zDiv
            }// end of i>=0 && j>=0 && k>=0 
         }
      }
   }// end of i,j,k

   if(type == 1){
      writeLogNum("wallmesh.log","MAXIMUM FACES IN BD BOX (stator)", maxFaces);
   }
   else{
      writeLogNum("wallmesh.log","MAXIMUM FACES IN BD BOX (rotor)", maxRotFaces);
   }
}

/*
* Go through all bdBox cell array and Insert faces into sorted face 
* array according to ascending order of cell index
*/
void insertToSortedFace(int xDv, int yDv, int zDv,
   int *con1, int *con2, int *con3, int *nCon1, int *nCon2, int *nCon3, int *cPtr, int type){
   int count = 0;
   for(int k=0; k<zDv; k++){
      for(int j=0; j<yDv; j++){
         for(int i=0; i<xDv; i++){
            int cI = i + j*xDv + k*xDv*yDv;
            // if(bdBox[cI].totalFacesR > 0){
            //    writeLogNum("wallmesh.log"," insertToSortedFace - no of faces ",bdBox[cI].totalFacesR);
            // }
            cPtr[cI] = count;
            if(type == 1){
               for(int n=0; n<bdBox[cI].totalFaces; n++){
                  nCon1[count] = con1[bdBox[cI].face[n]];
                  nCon2[count] = con2[bdBox[cI].face[n]];
                  nCon3[count] = con3[bdBox[cI].face[n]];
                  count++;
               }
            }
            // else{
            //    for(int n=0; n<bdBox[cI].totalFacesR; n++){
            //       nCon1[count] = con1[bdBox[cI].faceR[n]];
            //       nCon2[count] = con2[bdBox[cI].faceR[n]];
            //       nCon3[count] = con3[bdBox[cI].faceR[n]];
            //       count++;
            //    }
            // }
         }
      }
   }
   free(con1);
   free(con2);
   free(con3);
   for(int k=0; k<zDv; k++){
      for(int j=0; j<yDv; j++){
         for(int i=0; i<xDv; i++){
            int cI = i + j*xDv + k*xDv*yDv;
            if(type == 1){
               free(bdBox[cI].face);
            }
            // else{
            //    free(bdBox[cI].faceR);
            // }
         }
      }
   }
}


/*
* Calculate minimum and maximum face area
*/
void getMinMaxFaceArea(int noOfF, double *vX, double *vY, double *vZ,
      int *con1, int *con2, int *con3, int type){

   for(int i=0; i<noOfF; i++){
      double *node1 = allocateDoubleArray(DIM);
      double *node2 = allocateDoubleArray(DIM);
      double *node3 = allocateDoubleArray(DIM);
      getFace(node1,node2,node3, con1[i],con2[i],con3[i],vX,vY,vZ);
      double *n1n2 = allocateDoubleArray(DIM);
		double *n1n3 = allocateDoubleArray(DIM);
		vecSub(node1,node2,n1n2);
		vecSub(node1,node3,n1n3);
		double faceArea = triArea(n1n2,n1n3);
		minFaceArea = fmin(minFaceArea, faceArea);
		maxFaceArea = fmax(maxFaceArea, faceArea);
      minMaxEdgeLength(node1, node2, node3, type);
		free(n1n2);
		free(n1n3);
      free(node1);
      free(node2);
      free(node3);
   }
}

void minMaxEdgeLength(double *n1, double *n2, double *n3, int type){
    // Find min and mx of edge length
    double *tempLength = allocateDoubleArray(DIM);
    vecSub(n1,n2,tempLength);
    double tLength1 = vecMag(tempLength);
    vecSub(n2,n3,tempLength);
    double tLength2 = vecMag(tempLength);
    vecSub(n3,n1,tempLength);
    double tLength3 = vecMag(tempLength);
    if(type == 1){
    minEdgeLength = fmin(minEdgeLength,tLength1);
    maxEdgeLength = fmax(maxEdgeLength,tLength1);
    minEdgeLength = fmin(minEdgeLength,tLength2);
    maxEdgeLength = fmax(maxEdgeLength,tLength2);
    minEdgeLength = fmin(minEdgeLength,tLength3);
    maxEdgeLength = fmax(maxEdgeLength,tLength3);
    }

    else if(type == 2){
        minRotEdgeLength = fmin(minRotEdgeLength,tLength1);
        maxRotEdgeLength = fmax(maxRotEdgeLength,tLength1);
        minRotEdgeLength = fmin(minRotEdgeLength,tLength2);
        maxRotEdgeLength = fmax(maxRotEdgeLength,tLength2);
        minRotEdgeLength = fmin(minRotEdgeLength,tLength3);
        maxRotEdgeLength = fmax(maxRotEdgeLength,tLength3);
    }       

    free(tempLength); 
}

void getFace(double *n1, double *n2, double *n3, 
   int c1, int c2, int c3, double *vX, double *vY, double *vZ){
   n1[0] = vX[c1];
   n1[1] = vY[c1];
   n1[2] = vZ[c1];
   n2[0] = vX[c2];
   n2[1] = vY[c2];
   n2[2] = vZ[c2];
   n3[0] = vX[c3];
   n3[1] = vY[c3];
   n3[2] = vZ[c3];
}
