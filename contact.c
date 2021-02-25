#include "common.h"

/*Calculate overlap with contact surface
param:
parPos - particle center
dia - particle diameter
node1, node2, node3 - contact surface nodes
uVec - unit vector normal to surface
*/
double getOverlap(double *parPos, double dia, double *n1, double *n2, double *n3, double *uVec)
{
    double *v1 = allocateDoubleArray(DIM);     //vector running from node1 to particles center
    double *v2 = allocateDoubleArray(DIM);     //vector running from node1 to particles center
    double *ppDash = allocateDoubleArray(DIM); //vector running from particle center to projection
    vecSub(parPos, n1, v1);
    projVec(v1, uVec, v2, 0);
    vecSub(v1, v2, ppDash);

    double overlap = vecMag(ppDash) - 0.5 * dia;

    free(v1);
    free(v2);
    free(ppDash);
    return overlap;
}

bool inContact(int ip, int cnt)
{
    for(int i=0; i<demPart[ip].noOfCnt; i++)
    {
        if(cnt == demPart[ip].contList[i]){
            return true;
        }

    }
    return false;
}

void findLifterContact(int ip, double omega){
	int noOfLifters = 8;
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0;
    nFace = 2; //lifter contact
	double x0 = cylR - 0.5*lifterHeight;
	//double y0 = 0.0;
	
	for(int i=0; i<noOfLifters; i++){
		double rAng = i*2.0*M_PI/noOfLifters + rotorAngPosition;
		// Particle position in rotated coordinate system
		double xx = demPart[ip].pos[0] * cos(rAng) + demPart[ip].pos[2] * sin(rAng);
		double yy = demPart[ip].pos[1];
        double zz = demPart[ip].pos[2] * cos(rAng) - demPart[ip].pos[0] * sin(rAng);

    //         double cntPntX = cntPnt[0] * cos(-rAng) + cntPnt[2] * sin(-rAng);
    // double cntPntY = cntPnt[1];
    // double cntPntZ = cntPnt[2] * cos(-rAng) - cntPnt[0] * sin(-rAng);

        rVec[0] = -demPart[ip].pos[0];
        rVec[1] = 0.0;
        rVec[2] = -demPart[ip].pos[2];
        unitVec(rVec, uVec);
		// contact with top plane
		if(fabs(zz) <= 0.5*lifterWidth)
		{
            cntPnt[0] = x0 - 0.5*lifterHeight;
            cntPnt[1] = yy;
            cntPnt[2] = zz;
			contactPointInXYCoord(cntPnt, ip, rAng);		
		}
		else if(xx >= x0 - 0.5*lifterHeight){
			if(fabs(yy) - 0.5*demPart[ip].dia < 0.5*lifterWidth){
            cntPnt[0] = xx;
            cntPnt[1] = yy;
            cntPnt[2] = 0.5*lifterWidth*fabs(zz)/zz;
				contactPointInXYCoord(cntPnt, ip, rAng);
			}
		}
		else{
			// edge contact
            cntPnt[0] = x0 - 0.5*lifterHeight;
            cntPnt[1] = yy;
            cntPnt[2] = 0.5*lifterWidth*fabs(zz)/zz;
			contactPointInXYCoord(cntPnt, ip, rAng);
		}
	}
}

void findDrumContact(int ip, double omega)
{
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0;
    nFace = 1; //drum contact

    double radialDist = sqrt((demPart[ip].pos[0] * demPart[ip].pos[0]) + 
        (demPart[ip].pos[2] * demPart[ip].pos[2]));

    if(radialDist > (cylR - demPart[ip].dia*0.5)){
        rVec[0] = -demPart[ip].pos[0];
        rVec[1] = 0.0;
        rVec[2] = -demPart[ip].pos[2];
        unitVec(rVec, uVec);               
        cntPnt[0] = -uVec[0] * cylR;
        cntPnt[1] = demPart[ip].pos[1];
        cntPnt[2] = -uVec[2] * cylR;
        //printf("cntPnt %lf\n",cntPnt[2]/lengthFactor);
        contactPointInXYCoord(cntPnt, ip, 0.0);
    }   
}


void findEndWallContact(int ip, double omega)
{
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0;

    if(demPart[ip].pos[1] <= 0.5*demPart[ip].dia)
    {
        nFace = 3; //first end wall contact
        rVec[0] = 0.0;
        rVec[1] = demPart[ip].pos[1];
        rVec[2] = 0.0;
        unitVec(rVec, uVec);               
        cntPnt[0] = demPart[ip].pos[0];
        cntPnt[1] = 0.0;
        cntPnt[2] = demPart[ip].pos[2];
        //printf("cntPnt %lf\n",cntPnt[2]/lengthFactor);
        contactPointInXYCoord(cntPnt, ip, 0.0);        
    }
    else if(demPart[ip].pos[1] >= cylL - 0.5*demPart[ip].dia)
    {
        nFace = 4; //second end wall contact
        rVec[0] = 0.0;
        rVec[1] = -demPart[ip].pos[1];
        rVec[2] = 0.0;
        unitVec(rVec, uVec);               
        cntPnt[0] = demPart[ip].pos[0];
        cntPnt[1] = cylL;
        cntPnt[2] = demPart[ip].pos[2];
        //printf("cntPnt %lf\n",cntPnt[2]/lengthFactor);
        contactPointInXYCoord(cntPnt, ip, 0.0);        
    }
}

void neighbourContactForce(int pI)
{
    int end;
    if (pI == np-1)
    {
        end = nListSize;
    }
    else
    {
        end = nPtr[pI + 1];
    }
    for (int n = nPtr[pI]; n < end; n++)
    {
        int jp = neighList[n];
        double gap = getCenterDist(pI, jp) - (demPart[pI].dia + demPart[jp].dia) * 0.5;
        vecSub(demPart[pI].pos, demPart[jp].pos, ijVec);
        unitVec(ijVec,uVec);
        //ppElectForce(pI, jp, gap, uVec);
        if(gap < 100.e-9*lengthFactor)//activate vanderwaal force when gap<100nm
        {
            //ppVWForce(pI, jp, gap, uVec);
            //ppCapillaryForce(pI, jp, gap, uVec);
            if (gap < 0.0)
            {
                partContactForce(pI, jp, -gap, uVec);
                if(contactExist(pI, jp) == 0) //not in contact
                {
                    impactVelocityPP(pI, jp, uVec);
                    //ppCharge(pI, jp); //always place this line below impactVelocityPP
                    addPPContact(pI, demPart[pI].impAngPre[demPart[pI].noOfCnt], 
                        demPart[pI].impVelPre[demPart[pI].noOfCnt], jp);
                    // if jp is in contact with wall: ip wall contact = jp wall contact
                    if(demPart[jp].wallCnt != 0){
                        demPart[pI].wallCnt = demPart[jp].wallCnt;
                    }           
                }
                //if(!inContact(pI, jp)){
                    //double impAng = 0.0;
                    // vecAdd(demPart[pI].vel, demPart[jp].vel, tempVel);
                    // double impVel = vecMag(tempVel);
                    //addPPContact(pI, demPart[pI].impAngPre, demPart[pI].impVelPre, jp);             
                    
                //}
            }

        }
        if(gap > contactGap*demPart[pI].dia){
            deleteContact(pI, jp);
            deleteContact(jp, pI);
        }
    }
}


/*void findContactWallContact(Particle *p)
{
    int n, nearNode, noOfF;
    noOfF = 0;
    int *faceArr = allocateIntArray(10);
    double *nearNDPos = allocateDoubleArray(DIM);
    face_t f;
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    Thread *tf;
    double dist = 1.0e6;
    c_node_loop(c, tc, n)
    {        
        Node *v = C_NODE(c, tc, n);
        double *tempPos = allocateDoubleArray(DIM);
        double *tempVec = allocateDoubleArray(DIM);
        tempPos[0] = NODE_X(v);
        tempPos[1] = NODE_Y(v);
        tempPos[2] = NODE_Z(v);

        vecSub(demPart[p->part_id].pos, tempPos, tempVec);
        if(vecMag(tempVec) < dist)
        {
            noOfF = 0;
            nearNDPos[0] = NODE_X(v);
            nearNDPos[1] = NODE_Y(v);
            nearNDPos[2] = NODE_Z(v);
            nearNode = N_UDMI(v, 0);
            for(int i=0; i<N_UDMI(v, 1); i++)
            {
                int fid = N_UDMI(v, 2+i);
                faceArr[noOfF] = fid;
                noOfF++;
            }
            //v2 = v;
            dist = vecMag(tempVec);
        }
        free(tempVec);
        free(tempPos);
    
    }

    for(int i=0; i<noOfF; i++) //for number of faces
    {
        int nd1 = fcArray[faceArr[i]].wNode[0];
        int nd2 = fcArray[faceArr[i]].wNode[1];
        int nd3 = fcArray[faceArr[i]].wNode[2];

        if(!fcArray[faceArr[i]].rot)
        {
            //checkFaceInTriangle(p, ndArray[nd1].coord, ndArray[nd2].coord, ndArray[nd3].coord, 0.0, 0.0);
        }
    } 
    free(nearNDPos);
    free(faceArr);
}*/

void checkFaceInTriangle(int ip, double *node1, double *node2, double *node3, double rAng, double omega)
{
    rotOmega[0] = 0.0;
    rotOmega[1] = omega;
    rotOmega[2] = 0.0; 

    //int ip = p->part_id;
    double *pDash = allocateDoubleArray(DIM);
    pDash[0] = demPart[ip].pos[axis1] * cos(rAng) + demPart[ip].pos[axis2] * sin(rAng);
    pDash[1] = demPart[ip].pos[1];
    pDash[2] = demPart[ip].pos[axis2] * cos(rAng) - demPart[ip].pos[axis1] * sin(rAng);
        //--------- Face in triangle ----
    double *dVec = allocateDoubleArray(DIM);

        if (getProjectOnFace(pDash, node1, node2, node3, dVec))
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
    free(pDash);

}

void contactPointInXYCoord(double *cntPnt, int ip, double rAng){
    //contact point in XY coordinate system
    double cntPntX = cntPnt[0] * cos(-rAng) + cntPnt[2] * sin(-rAng);
    double cntPntY = cntPnt[1];
    double cntPntZ = cntPnt[2] * cos(-rAng) - cntPnt[0] * sin(-rAng);
    torqueRadius = sqrt(cntPntX*cntPntX + cntPntZ*cntPntZ);
 
    jpRVec[0] = -cntPntX;
    jpRVec[1] = -cntPntY;
    jpRVec[2] = -cntPntZ;
    double *dVecXY = allocateDoubleArray(DIM);
    dVecXY[0] = demPart[ip].pos[0] - cntPntX;
    dVecXY[1] = demPart[ip].pos[1] - cntPntY;
    dVecXY[2] = demPart[ip].pos[2] - cntPntZ;

    unitVec(dVecXY, uVec);
 
    findContactFromMesh(ip, dVecXY);
    free(dVecXY);
}

/*Find particle-wall contact and calculate force using triangular boundary Mesh*/
void findContactFromMesh(int ip, double *dVec)
{
    double gap = vecMag(dVec) - 0.5 * demPart[ip].dia;
    // printf("GAP %lf\n",gap/lengthFactor);
    // exit(0);
    demPart[ip].minWallDist = fmin(gap, demPart[ip].minWallDist);
    //pwElectForce(ip, gap, uVec);
    
        if (gap < 100.e-9 * lengthFactor)
        {
            pWallVWForce(ip, gap, uVec);
            //pwCapillaryForce(ip, gap, uVec);
            if (gap < 0) //If contact exists calculate contact force
            {
                //printf("GAP %lf\n",gap/lengthFactor);
                demPart[ip].wallCnt = nFace;
                surfaceContactForce(ip, -gap, uVec);
                if(demPart[ip].impactExist == false){
                    impactVelocity(ip, uVec);
                    //demPart[ip].trapped = true;
                    //pwCharge(ip);
                    demPart[ip].impactExist = true;      
                    addImpact(ip, demPart[ip].preRelAngWall, demPart[ip].preRelVelWall, nFace);
                }
            }
        }
    
}

/*void checkXContact(Particle *p, double xMin, double xMax)
{
    //Contact with xMin
    double gap = demPart[p->part_id].pos[0] - xMin - demPart[p->part_id].dia * 0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;

    if (gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100.e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }
    //Contact with xMax
    gap = xMax - demPart[p->part_id].pos[0] - demPart[p->part_id].dia * 0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;
    if (gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100.e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }
    //if(gap < demPart[p->part_id].dia*0.05){pWallVWForce(p->part_id, gap, uVec);}
}*/

/*
void checkZContact(Particle *p, double zMin, double zMax)
{
    //Contact with zMin
    double gap = -(zMin - demPart[p->part_id].pos[2]) - demPart[p->part_id].dia * 0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = 1.0;
    if (gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }
    //Contact with z=88mm wall
    gap = zMax - demPart[p->part_id].pos[2] - demPart[p->part_id].dia * 0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;
    if (gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100.e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }
}

void checkYContact(Particle *p, double yMin, double yMax)
{
    // Contact with yMin
    double gap = -(yMin - demPart[p->part_id].pos[1]) - demPart[p->part_id].dia * 0.5;
    uVec[0] = 0.0;
    uVec[1] = 1.0;
    uVec[2] = 0.0;

    if (gap < 0) //If contact exists calculate contact force
    {
        //demPart[p->part_id].force[1] += -demPart[p->part_id].mass;
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100.e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }

    // Contact with yMax
    gap = yMax - (demPart[p->part_id].pos[1] + demPart[p->part_id].dia * 0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;

    if (gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }
    if (gap < 100.e-9 * lengthFactor)
    {
        pWallVWForce(p->part_id, gap, uVec);
    }
}*/
