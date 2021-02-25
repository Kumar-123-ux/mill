#include "common.h"

/* Return particle volume*/
double partVol(int p)
{
    return (4.0 / 3.0) * PI * pow(demPart[p].dia * 0.5, 3);
}

/*Find solid fraction within cell radius*/
double solidFraction(int ip)
{
    //Find cell center
    int iIndex = floor((demPart[ip].pos[0] - xmin) / domainDx);
    int jIndex = floor((demPart[ip].pos[1] - ymin) / domainDy);
    int kIndex = floor((demPart[ip].pos[2] - zmin) / domainDz);
    int cellIndex = iIndex + jIndex * xDiv + kIndex * xDiv * yDiv;

    double cX = iIndex * domainDx + 0.5 * domainDx;
    double cY = jIndex * domainDy + 0.5 * domainDy;
    double cZ = kIndex * domainDz + 0.5 * domainDz;
    double vol = 0.0;
    double totVol = 0.0;
    for (int i = 0; i < bdBox[cellIndex].noOfParticles; i++)
    {
        int jp = bdBox[cellIndex].parts[i];
        double r = 0.5 * demPart[jp].dia;
        double jpX = demPart[jp].pos[0];
        double jpY = demPart[jp].pos[1];
        double jpZ = demPart[jp].pos[2];

        double rR = sqrt((jpX - cX) * (jpX - cX) + (jpY - cY) * (jpY - cY) + (jpZ - cZ) * (jpZ - cZ));

        if (rR <= (0.5 * domainDx - r))
        {
            vol = 0.8 * partVol(jp);
        }
        else if (rR > (cellRadius - r))
        {
            vol = partVol(jp);
        }
        else if (rR <= (cellRadius - r) && rR > (domainDx * 0.5 - r))
        {
            vol = 0.5 * partVol(jp);
        }
        totVol += vol;
    }
    //writeLogNum("logfile3.log","VOL ",totVol/(lengthFactor*lengthFactor*lengthFactor));
    //bdBox[cellIndex].fluidVolF = totVol/((4.0/3.0)*PI*pow(cellRadius,3));
    return totVol / ((4.0 / 3.0) * PI * pow(cellRadius, 3));
}

/*
  Search through all the faces in neighbour cells and return the nearest face
*/
void getNearFace(int ip, int *con1, int *con2, int *con3,
                 double *vX, double *vY, double *vZ, int *cPtr,
                 double dDx, double dDy, double dDz, int xD, int yD, int zD, 
                 int noOfF, double rAng, double omega)
{
    //Find cell center
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0; 

    int minCI = 1e5;
    int maxCI = -100;
    
    double *pDash = allocateDoubleArray(DIM);
    pDash[0] = demPart[ip].pos[axis1] * cos(rAng) + demPart[ip].pos[axis2] * sin(rAng);
    pDash[1] = demPart[ip].pos[1];
    pDash[2] = demPart[ip].pos[axis2] * cos(rAng) - demPart[ip].pos[axis1] * sin(rAng);

    //if(iter%2 == 0){
    int iIndex = floor((pDash[0] - xmin) / dDx);
    int jIndex = floor((pDash[1] - ymin) / dDy);
    int kIndex = floor((pDash[2] - zmin) / dDz);
    double prevDist = 1.e3;
    double dist = 0.;
    
    int end;
    int cI = iIndex + jIndex * xD + kIndex * xD * yD;
    minCI = fmin(cI, minCI);
    maxCI = fmax(cI, maxCI);

    if ((cI + 1) > xD * yD * zD)
    {
        end = noOfF;
    }
    else
    {
        end = cPtr[cI + 1];
    }

    for (int n = cPtr[cI]; n < end; n++)
    {
        double *nd1 = allocateDoubleArray(DIM);
        double *nd2 = allocateDoubleArray(DIM);
        double *nd3 = allocateDoubleArray(DIM);
        nd1[0] = vX[con1[n]];
        nd1[1] = vY[con1[n]];
        nd1[2] = vZ[con1[n]];
        nd2[0] = vX[con2[n]];
        nd2[1] = vY[con2[n]];
        nd2[2] = vZ[con2[n]];
        nd3[0] = vX[con3[n]];
        nd3[1] = vY[con3[n]];
        nd3[2] = vZ[con3[n]];

        //--------- Face in triangle ----
        double *dVec = allocateDoubleArray(DIM);

        if (getProjectOnFace(pDash, nd1, nd2, nd3, dVec))
        {
            nFace = 1;
           //double *uVec = allocateDoubleArray(DIM);
            double *cntPnt = allocateDoubleArray(DIM);
            //Since dVec is point towards center we take minus to get contact point
            double cntPntX = -dVec[0] + pDash[0];
            double cntPntY = -dVec[1] + pDash[1];
            double cntPntZ = -dVec[2] + pDash[2];
            
            cntPnt[0] = cntPntX * cos(-rAng) + cntPntZ * sin(-rAng);
            cntPnt[1] = cntPntY;
            cntPnt[2] = cntPntZ * cos(-rAng) - cntPntX * sin(-rAng);

            jpRVec[0] = -cntPnt[0];
            jpRVec[1] = -cntPnt[1];
            jpRVec[2] = -cntPnt[2];

            double *dVecXY = allocateDoubleArray(DIM);
            dVecXY[0] = demPart[ip].pos[0] - cntPnt[0];
            dVecXY[1] = demPart[ip].pos[1] - cntPnt[1];
            dVecXY[2] = demPart[ip].pos[2] - cntPnt[2];
            unitVec(dVecXY, uVec);
            findContactFromMesh(ip, dVecXY);

            free(dVecXY);
            free(cntPnt);
        }
    
        free(dVec);
        free(nd1);
        free(nd2);
        free(nd3);
        
        //---------end of face in triangle
    }

    free(pDash);
}

void outerWall(int ip){
    
    double cylR = 9.5 * 1e-3 * lengthFactor;
    double endCap1 = -3.7 * 1e-3 * lengthFactor;
    double endCap2 = 3.7 * 1e-3 * lengthFactor;

    double radialDist = sqrt((demPart[ip].pos[0] * demPart[ip].pos[0]) + 
        (demPart[ip].pos[2] * demPart[ip].pos[2]));

    if(radialDist > (cylR - demPart[ip].dia*0.5)){
        rVec[0] = -demPart[ip].pos[0];
        rVec[1] = 0.0;
        rVec[2] = -demPart[ip].pos[2];
        unitVec(rVec, uVec);               
        cntPnt[0] = uVec[1] * cylR;
        cntPnt[1] = 0.0;
        cntPnt[2] = uVec[2] * cylR;
        contactPointInXYCoord(cntPnt, ip, 0.0);
    }
    if(demPart[ip].pos[1] < (endCap1 + demPart[ip].dia)) // bottom end cap
    {
        rVec[0] = 0.0;
        rVec[1] = -demPart[ip].pos[1];
        rVec[2] = 0.0;
        unitVec(rVec, uVec);               
        cntPnt[0] = 0.0;
        cntPnt[1] = uVec[1] * endCap1;
        cntPnt[2] = 0.0;
        contactPointInXYCoord(cntPnt, ip, 0.0);
    }
    if(demPart[ip].pos[1] > (endCap2 - demPart[ip].dia*0.5)) //top end cap
    {
        rVec[0] = 0.0;
        rVec[1] = -demPart[ip].pos[1];
        rVec[2] = 0.0;
        unitVec(rVec, uVec);               
        cntPnt[0] = 0.0;
        cntPnt[1] = uVec[1] * endCap2;
        cntPnt[2] = 0.0;
        contactPointInXYCoord(cntPnt, ip, 0.0);
    }
}

/*void capsuleFace(int ip, double rAng, double omega)
{
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0; 

    double pX = 0.5*cylL+cylR;
    double pZ = 0.0;
    //virtual sphear at center of hole 1
    holeCenter1[0] = pX*cos(-rAng) + pZ*sin(-rAng);
    holeCenter1[2] = pZ*cos(-rAng) - pX* sin(-rAng);

    pX = -(0.5*cylL+cylR);
    pZ = 0.0;
    //virtual sphear at center of hole 2
    holeCenter2[0] = pX*cos(-rAng) + pZ*sin(-rAng);
    holeCenter2[2] = pZ*cos(-rAng) - pX* sin(-rAng);

    double *pDash = allocateDoubleArray(DIM); //particle position in rotated coordinate
    pDash[0] = demPart[ip].pos[0] * cos(rAng) + demPart[ip].pos[2] * sin(rAng);
    pDash[1] = demPart[ip].pos[1];
    pDash[2] = demPart[ip].pos[2] * cos(rAng) - demPart[ip].pos[0] * sin(rAng);

    //double cylR = 2.0 * 1e-3 * lengthFactor;
    //double cylL = 10.0 * 1e-3 * lengthFactor;
    //double holeR = 0.42 * 1e-3 * lengthFactor;
    double edgeDia = 0.140*1e-3 * lengthFactor;

    double dist1 = sqrt(pow((pDash[0]-holeCenter1[0]),2)
        +pow((pDash[1]-holeCenter1[1]),2)
        +pow((pDash[2]-holeCenter1[2]),2));
    double dist2 = sqrt(pow((pDash[0]-holeCenter2[0]),2)
        +pow((pDash[1]-holeCenter2[1]),2)
        +pow((pDash[2]-holeCenter2[2]),2));    
    if(dist1 < holeR || dist2 < holeR){
        holePart[ip] = 1;
    }
    // chack for contact with cylinder
    if (pDash[0] <= cylL * 0.5 && pDash[0] >= -cylL * 0.5)
    {
        nFace = -1;
        rVec[0] = 0.0;
        rVec[1] = pDash[1];
        rVec[2] = pDash[2];

        unitVec(rVec, uVec);
        //cntPnt = allocateDoubleArray(DIM);
        //double cntPnt[DIM];
        cntPnt[0] = pDash[0];
        cntPnt[1] = uVec[1] * cylR;
        cntPnt[2] = uVec[2] * cylR;

        contactPointInXYCoord(cntPnt, ip, rAng);
    }
    else if(pDash[0] > cylL * 0.5 || pDash[0] < -cylL * 0.5) //contact with end cap
    {
        // check for contact with end cap
        
        double radDist = sqrt(pDash[1] * pDash[1] + pDash[2] * pDash[2]);
        if(radDist > holeR)
        {
            //Can have a contact with cap ends or cap hole edge
            if (pDash[0] > cylL * 0.5)
            {
                nFace = -2;
                capCent[0] = cylL * 0.5;
                capCent[1] = 0.0;
                capCent[2] = 0.0;
            }
            //Can have a contact with cap ends or cap hole edge
            else if (pDash[0] < -cylL * 0.5)
            {
                nFace = -3;
                capCent[0] = -cylL * 0.5;
                capCent[1] = 0.0;
                capCent[2] = 0.0;
            }
            vecSub(pDash, capCent, rVec); //vector from cap center to particle
            unitVec(rVec, uVec);
            sclVecMult(cylR, uVec, rVec); //vector from cap center to contact point
            vecAdd(capCent, rVec, cntPnt);
            contactPointInXYCoord(cntPnt, ip, rAng);
        }
        if(radDist > 0.5 * demPart[ip].dia) //inside the cap hole
        {
            tempUV[0] = 0.0;
            tempUV[1] = pDash[1];
            tempUV[2] = pDash[2];
            unitVec(tempUV, tempUV);

            if (pDash[0] > cylL * 0.5)
            {
                nFace = -4;
                edgeCent[0] = 0.5 * cylL + cylR + 0.5 * edgeDia;
                edgeCent[1] = tempUV[1] * holeR;
                edgeCent[2] = tempUV[2] * holeR;
                // double dist = sqrt(pow((pDash[0] - edgeCent[0]),2)+
                //     pow((pDash[1] - edgeCent[1]),2)+
                //     pow((pDash[2] - edgeCent[2]),2));
                // if(dist < holeR){
                //     holePart[ip] = 1;
                // }
            }
            else if (pDash[0] < -cylL * 0.5)
            {
                nFace = -5;
                edgeCent[0] = -(0.5 * cylL + cylR + 0.5 * edgeDia);
                edgeCent[1] = tempUV[1] * holeR;
                edgeCent[2] = tempUV[2] * holeR;
                // double dist = sqrt(pow((pDash[0] - edgeCent[0]),2)+
                //     pow((pDash[1] - edgeCent[1]),2)+
                //     pow((pDash[2] - edgeCent[2]),2));
                // if(dist < holeR){
                //     holePart[ip] = 1;
                // }      
            }

            vecSub(pDash, edgeCent, rVec);
            unitVec(rVec, uVec);
            sclVecMult(0.5*edgeDia, uVec, rVec);
            vecAdd(edgeCent, rVec, cntPntEdge);
            contactPointInXYCoord(cntPntEdge, ip, rAng);
        }
    }

    free(pDash);
}*/

double triArea(double *v1, double *v2)
{
    double a = 0;
    //Cross product
    double *v1v2 = allocateDoubleArray(DIM);
    crossProd(v1, v2, v1v2);
    a = 0.5 * vecMag(v1v2);
    free(v1v2);
    return a;
}

/* projection to surface*/
void getProjection(double *p, double *n1, double *n2, double *n3, double *dVec)
{
    //writeLog3Num("wallmesh.log"," particle ",p[0],p[1],p[2]);
    double *n1n2 = allocateDoubleArray(DIM);
    double *n1n3 = allocateDoubleArray(DIM);
    double *n1p = allocateDoubleArray(DIM);
    vecSub(n2, n1, n1n2);
    vecSub(n3, n1, n1n3);
    vecSub(p, n1, n1p);

    //Temp unit vector
    double *uVec = allocateDoubleArray(DIM);
    getUnitVector(n1, n2, n3, uVec);
    double *pj = allocateDoubleArray(DIM);
    projVec(n1p, uVec, pj, 0);

    //Projection to particle vector
    vecSub(pj, n1p, dVec);
    sclVecMult(-1.0, dVec, dVec);
    free(uVec);
    free(pj);
    free(n1n2);
    free(n1n3);
    free(n1p);
}

/* If projection is inside the triangle return 1 else return 0*/
bool getProjectOnFace(double *p, double *n1, double *n2, double *n3, double *dVec)
{
    //writeLog3Num("wallmesh.log"," particle ",p[0],p[1],p[2]);
    double *n1n2 = allocateDoubleArray(DIM);
    double *n2n3 = allocateDoubleArray(DIM);
    double *n3n1 = allocateDoubleArray(DIM);
    double *n1p = allocateDoubleArray(DIM);
    vecSub(n2, n1, n1n2);
    vecSub(n3, n2, n2n3);
    vecSub(n1, n3, n3n1);
    vecSub(p, n1, n1p);

    //Temp unit vector
    double *uVec = allocateDoubleArray(DIM);
    getUnitVector(n1, n2, n3, uVec);
    double *pj = allocateDoubleArray(DIM);
    projVec(n1p, uVec, pj, 0);

    //Projection to particle vector
    vecSub(n1p, pj, dVec); //dVec pointing to particle center
    //sclVecMult(-1.0, dVec, dVec);

    //vector from projection to n2 node
    double *n2Proj = allocateDoubleArray(DIM);
    vecSub(pj,n1n2, n2Proj);
    //vector from projection to n3 node
    double *n3Proj = allocateDoubleArray(DIM);
    vecAdd(pj, n3n1, n3Proj);
    //writeLog3Num("wallmesh.log","proj face ",pj[0],pj[1],pj[2]);

    double a1 = triArea(n1n2, pj);
    double a2 = triArea(n2n3, n2Proj);
    double a3 = triArea(n3n1, n3Proj);

    double *n1n3 = allocateDoubleArray(DIM);
    sclVecMult(-1., n3n1, n1n3);
    double totalArea = triArea(n1n2, n1n3);

    free(n2Proj);
    free(n3Proj);
    free(n1p);
    free(pj);
    free(uVec);
    free(n1n2);
    free(n2n3);
    free(n3n1);
    free(n1n3);

    //double tolerence = minFaceArea * (0.1 / 100.) * pow(lengthFactor, 2);
    //if(fabs(a1+a2+a3 - totalArea) < totalArea*1.e-3)
    if(fabs(a1+a2+a3 - totalArea) < totalArea*0.001)
    //if (fabs(a1 + a2 + a3 - totalArea) < tolerence)
    {
        return true;
    }
    return false;
}

/* Return center distance of two particles
param:
p1 - particle 1
p2 - particle 2 */
double getCenterDist(int ip, int jp)
{
    double ipX = demPart[ip].pos[0];
    double ipY = demPart[ip].pos[1];
    double ipZ = demPart[ip].pos[2];
    double jpX = demPart[jp].pos[0];
    double jpY = demPart[jp].pos[1];
    double jpZ = demPart[jp].pos[2];

    double val = sqrt((ipX - jpX) * (ipX - jpX) + (ipY - jpY) * (ipY - jpY) + (ipZ - jpZ) * (ipZ - jpZ));
    return val;
}

double getDist(double *v1, double *v2){
    double dist = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]));
    return dist;
}

/* Add two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecAdd(double *v1, double *v2, double *vec)
{
    vec[0] = v1[0] + v2[0];
    vec[1] = v1[1] + v2[1];
    vec[2] = v1[2] + v2[2];
}

/* Substract two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecSub(double *v1, double *v2, double *vec)
{
    vec[0] = v1[0] - v2[0];
    vec[1] = v1[1] - v2[1];
    vec[2] = v1[2] - v2[2];
}

/* Find relative velocity magnitude between two particles */
double relVel(int ip, int jp)
{
    double *relVelVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel, demPart[jp].vel, relVelVec);
    double vel = vecMag(relVelVec);
    free(relVelVec);
    return vel;
}

/* Find vector cross product
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void crossProd(double *v1, double *v2, double *vec)
{
    double temp1 = v1[1] * v2[2];
    double temp2 = v2[1] * v1[2];
    vec[0] = temp1 - temp2;

    temp1 = v2[0] * v1[2];
    temp2 = v1[0] * v2[2];
    vec[1] = temp1 - temp2;

    temp1 = v1[0] * v2[1];
    temp2 = v2[0] * v1[1];
    vec[2] = temp1 - temp2;
}

/*
Dot product of two vectors
param:
v1 - vector 1
v2 - vector 2
*/
double dotProduct(double *v1, double *v2)
{
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

/*
Find normal unit vector to the surface defined by vector v1v2 and v2v3
*/
void getUnitVector(double *v1, double *v2, double *v3, double *uVec)
{
    double *v1v2 = allocateDoubleArray(DIM);
    double *v1v3 = allocateDoubleArray(DIM);

    vecSub(v1, v2, v1v2);
    vecSub(v1, v3, v1v3);
    crossProd(v1v2, v1v3, uVec);
    unitVec(uVec, uVec);

    free(v1v2);
    free(v1v3);
}

/* Find unit vector
param: 
v - input vector
vec - unit vector */
void unitVec(double *v, double *vec)
{
    double vMag = vecMag(v);
    if (vMag > 0.0)
    {
        double temp1 = v[0] / vMag;
        double temp2 = v[1] / vMag;
        double temp3 = v[2] / vMag;
        vec[0] = temp1;
        vec[1] = temp2;
        vec[2] = temp3;
    }
    else
    {
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    }
}

/* Find magnitude of a vector*/
double vecMag(double *vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

/* Multiply input vector by a scaler
param: 
scl - scaler
vec - vector to be multiplied*/
void sclMult(double scl, double *vec)
{
    vec[0] = scl * vec[0];
    vec[1] = scl * vec[1];
    vec[2] = scl * vec[2];
}

/* Multiply input vector by a scaler and assign to a vector
param: 
scl - scaler
inVec - input vector to be multiplied
outVec - reusltant output vector*/
void sclVecMult(double scl, double *inVec, double *outVec)
{
    outVec[0] = scl * inVec[0];
    outVec[1] = scl * inVec[1];
    outVec[2] = scl * inVec[2];
}

/* Returns the project vector on the plane defined by normal unit vector
param:
v - input vector
n - unit vector
vec - resultant project vector
type - either 0 or 1, 0 gives project vector, 1 gives project vector scaled by input vector */
void projVec(double *v, double *n, double *vec, int type)
{
    double *tV1, *tV2;
    tV1 = (double*)malloc(DIM * sizeof(double));
    tV2 = (double*)malloc(DIM * sizeof(double));
    if (type == 0)
    {
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        //printf("tV1 %lf,%lf,%lf\n ",tV1[0],tV1[1],tV1[2]);
        //printf("tV2 %lf,%lf,%lf\n ",tV2[0],tV2[1],tV2[2]);
        crossProd(tV2, tV1, vec);
        //printf("Type 0\n");
    }
    else
    {
        double *tV3 = (double*)malloc(DIM * sizeof(double));
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        crossProd(tV2, tV1, tV3);
        double temp;
        if ((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) != 0.0)
        {
            temp = sqrt((tV3[0] * tV3[0] + tV3[1] * tV3[1] + tV3[2] * tV3[2]) /
                        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
        }
        else
        {
            temp = 0.0;
        }
        vec[0] = tV3[0] * temp;
        vec[1] = tV3[1] * temp;
        vec[2] = tV3[2] * temp;
        free(tV3);
    }
    free(tV1);
    free(tV2);
}
