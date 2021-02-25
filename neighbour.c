#include "common.h"

/*
Check for existing neighbour
*/
// int insertable(int ip, int jp)
// {
//     for (int i = 0; i < demPart[ip].noOfNeigh; i++)
//     {
//         if (demPart[ip].neigh[i] == jp)
//         {
//             return 0;
//         }
//     }
//     return 1;
// }

/*
Add a new neighbour to the neighbour list
*/
// void addNeighbour(int ip, int jp)
// {
//     demPart[ip].neigh[demPart[ip].noOfNeigh] = jp;
//     demPart[ip].noOfNeigh++;
// }

/*
Add a new contact to the contact list
*/
// void addContact(int ip, int cnt)
// {
//     demPart[ip].contList[demPart[ip].noOfCnt] = cnt;
//     demPart[ip].noOfCnt++;
// }


/*
Delete all contacts except nFace
*/
// void deleteAllContacts(int ip, int nF)
// {
//     for (int i = 0; i < demPart[ip].noOfWallCnt; i++)
//     {
//         int wall = demPart[ip].wallCntList[i];
//         if (nF != wall)
//         {
//             demPart[ip].wallCntList[i] = demPart[ip].wallCntList[demPart[ip].noOfWallCnt - 1];
//             demPart[ip].noOfWallCnt--;
//         }
//     }
// }

/*
Add a new contact to the contact list
*/
void addPPContact(int ip, double impAng, double impVel, int jp)
{
    //ip particle
    impactPartPP1[ppTotalImpacts] = ip;
    impactPartPP2[ppTotalImpacts] = jp;
    ppImpactVel[ppTotalImpacts] = impVel/velocityFactor;
    ppImpactAng[ppTotalImpacts] = impAng;
    ppParDia1[ppTotalImpacts] = 1.e3*demPart[ip].dia/lengthFactor;
    ppParDia2[ppTotalImpacts] = 1.e3*demPart[jp].dia/lengthFactor;
    ppImpactPosX[ppTotalImpacts] = demPart[ip].pos[0]/lengthFactor;
    ppImpactPosY[ppTotalImpacts] = demPart[ip].pos[1]/lengthFactor;
    ppImpactPosZ[ppTotalImpacts] = demPart[ip].pos[2]/lengthFactor;
    ppTotalImpacts++;

    demPart[ip].contList[demPart[ip].noOfCnt] = jp;
    demPart[ip].noOfCnt++;
    demPart[ip].totalParColl++;

    //jp particle
    // impactPartPP[ppTotalImpacts] = jp;
    // ppImpactVel[ppTotalImpacts] = impVel/velocityFactor;
    // ppImpactAng[ppTotalImpacts] = impAng;
    // ppParDia[ppTotalImpacts] = 1.e3*demPart[jp].dia/lengthFactor;
    // ppImpactPosX[ppTotalImpacts] = demPart[jp].pos[0]/lengthFactor;
    // ppImpactPosY[ppTotalImpacts] = demPart[jp].pos[1]/lengthFactor;
    // ppImpactPosZ[ppTotalImpacts] = demPart[jp].pos[2]/lengthFactor;

    demPart[jp].contList[demPart[jp].noOfCnt] = ip;
    demPart[jp].noOfCnt++;
    demPart[jp].totalParColl++;

    if(ppTotalImpacts > MAXIMPACTS - 100){
        dumpPPImpact("pp-impact.dat");
    }

    if(demPart[ip].noOfCnt > CNTSIZE){
        //writeLogNum("logfile10.log"," demPart[ip].noOfCnt > CNTSIZE ", demPart[ip].noOfCnt);
        demPart[ip].noOfCnt = 0; //reset
    }
    if(demPart[jp].noOfCnt > CNTSIZE){
        //writeLogNum("logfile10.log"," demPart[jp].noOfCnt > CNTSIZE ", demPart[jp].noOfCnt);
        demPart[jp].noOfCnt = 0; //reset
    }
}

/*
* Add a new particle-wall impact
*/
void addImpact(int ip, double impAng, double impVel, int surface){
    demPart[ip].wallCnt = surface;
    impactPart[totalImpacts] = ip;
    impactSurface[totalImpacts] = surface;
    impactVel[totalImpacts] = impVel/velocityFactor;
    impactAng[totalImpacts] = impAng;
    parDia[totalImpacts] = 1.e3*demPart[ip].dia/lengthFactor;
    impactPosX[totalImpacts] = demPart[ip].pos[0]/lengthFactor;
    impactPosY[totalImpacts] = demPart[ip].pos[1]/lengthFactor;
    impactPosZ[totalImpacts] = demPart[ip].pos[2]/lengthFactor;
    totalImpacts++;

    if(totalImpacts > MAXIMPACTS - 100){
        //dumpImpact("wall-impact.dat");
    }
    // if(totalImpacts > MAXIMPACTS){
    //     impactVel = (int *)doubleloc(impactVel, (maxImpactSize + IMPOFFSET)*sizeof(int));
    //     maxImpactSize += IMPOFFSET;
    //     //dumpImpact("wall-impact.dat");
    // }
}

/*
Add a new contact to the contact list
*/
// void addWallContact(int ip, int cnt)
// {
//     demPart[ip].wallCntList[demPart[ip].noOfWallCnt] = cnt;
//     demPart[ip].noOfWallCnt++;

//     //writeLogNum("logfile5.log","ADD ", cnt);
//     if (demPart[ip].noOfWallCnt > WALLCNTSIZE)
//     {
//         writeLogNum("logfile4.log", "demPart[ip].noOfWallCnt > WALLCNTSIZE", demPart[ip].noOfWallCnt);
//     }
// }

/*
Delete contact
*/
// void deleteWallContact(int ip, int cnt)
// {
//     for (int i = 0; i < demPart[ip].noOfWallCnt; i++)
//     {
//         int wall = demPart[ip].wallCntList[i];
//         if (cnt == wall)
//         {
//             //writeLogNum("logfile5.log","DELETE ", cnt);
//             demPart[ip].wallCntList[i] = demPart[ip].wallCntList[demPart[ip].noOfWallCnt - 1];
//             demPart[ip].noOfWallCnt--;
//             //break;
//         }
//     }
// }

/* Check for exisiting contact
If exist return true otherwise return false
*/
int contactExist(int ip, int cnt)
{
    for (int i = 0; i < demPart[ip].noOfCnt; i++)
    {
        if (demPart[ip].contList[i] == cnt)
        {
            return 1;
            //break;
        }
    }
    return 0;
}

/*
Delete contact jp contact from ip and vise versa
*/
void deleteContact(int ip, int cnt)
{
    for (int i = 0; i < demPart[ip].noOfCnt; i++)
    {
        int jp = demPart[ip].contList[i];
        if (cnt == jp)
        {
            demPart[ip].contList[i] = demPart[ip].contList[demPart[ip].noOfCnt - 1];
            demPart[ip].noOfCnt--;
            break;
        }
    }
}

/*
Delete neighbour
*/
// void deleteNeighbour(int ip, int jp)
// {
//     for (int i = 0; i < demPart[ip].noOfNeigh; i++)
//     {
//         int neigh = demPart[ip].neigh[i];
//         if (neigh == jp)
//         {
//             demPart[ip].neigh[i] = demPart[ip].neigh[demPart[ip].noOfNeigh - 1];
//             demPart[ip].noOfNeigh--;
//             break;
//         }
//     }
// }

/*
Insert particle to BdBox
param:
pI - particle index
cI - cell index
*/
void insertToBdBox(int ip, int cI)
{
    bdBox[cI].parts[bdBox[cI].noOfParticles] = ip;
    bdBox[cI].noOfParticles++;
    demPart[ip].prevCellIndex = cI;
    if (bdBox[cI].noOfParticles > 0.8*bdBox[cI].maxSize)
    {
        //bdBox[cI].parts = resize_array(bdBox[cI].parts, bdBox[cI].maxSize + 5*NO_OF_PARTICLES_IN_BDCELL);
        bdBox[cI].maxSize = bdBox[cI].maxSize + 5*NO_OF_PARTICLES_IN_BDCELL;
        //writeLogNum("logfile10.log","expand BDBOX ",bdBox[cI].maxSize);       
    }
}

/*
Delete particle from bdBox
param:
pI - particle index
cI = cell index
*/
void deleteParticle(int p, int cI)
{
    for (int i = 0; i < bdBox[cI].noOfParticles; i++)
    {
        int np = bdBox[cI].parts[i];
        if (p == np)
        {
            bdBox[cI].parts[i] = bdBox[cI].parts[bdBox[cI].noOfParticles - 1];
            bdBox[cI].noOfParticles -= 1;
            break;
        }
    }
}

/*
bool inBondList(int ip, int jp){
    // for(int i=0; i<demPart[ip].bond; i++){
    //     if(demPart[ip].bondList[i] == jp)
    //     {
    //         return true;
    //     }
    // }
    return false;
}*/

/*
void addToBondList(int ip, int jp){
    // demPart[ip].bondList[demPart[ip].bond] = jp;
    // demPart[ip].bond += 1;
}
*/
/*
* find attached drug particles on to carrier
*/
/*
void findAttachedParticles(int ip){
    int end;
    if (ip == np)
    {
        end = nListSize;
    }
    else
    {
        end = nPtr[ip + 1];
    }
    for (int n = nPtr[ip]; n < end; n++)
    {
        int jp = neighList[n];
        vecSub(demPart[ip].vel, demPart[jp].vel, tempVel);
        if(vecMag(tempVel)/velocityFactor < 0.001){
            demPart[ip].bond += 1;
        }
        //findAttachedParticles(jp, origin);
    }  
}*/

/*
void findNext(int ip, int origin)
{
    int end;
    if (ip == np)
    {
        end = nListSize;
    }
    else
    {
        end = nPtr[ip + 1];
    }
    for (int n = nPtr[ip]; n < end; n++)
    {
        int jp = neighList[n];
        double gap = getCenterDist(ip, jp) - (demPart[ip].dia + demPart[jp].dia) * 0.5;
        if (gap < 100.e-9 * lengthFactor)
        {
            if(!inBondList(origin, jp))
            {
                //demPart[origin].bond += 1;
                addToBondList(origin, jp);
            }
            findNext(jp, origin);
        }
    }
}*/

/*
* Neighobulist update without neighbour cells
*/
// void updateNeighbourList()
// {
//     int count = 0;
//     for(int ip=0; ip<np; ip++){
//         nPtr[ip] = count;
//         for (int jp = ip + 1; jp < np; jp++){
            
//             vecSub(demPart[ip].pos,demPart[jp].pos,tempDist);
                
//             if(vecMag(tempDist) < 1.55*0.5*(demPart[ip].dia + demPart[jp].dia)){
//                 neighList[count] = jp;
//                 count++;
//             }
//         }
//         demPart[ip].displacement = 0.0;
//         demPart[ip].bond = 0;
//         //if(demPart[ip].dia < 5.e-6*lengthFactor){
//         //findNext(ip, ip); // find bond paritcles
//         //}
//     }
//     nListSize = count;
//     if(nListSize > NEIGH_LIST_ARRAY_SZ){
//         writeLogNum("logfile10.log", "nListSize > NEIGH_LIST_ARRAY_SZ ", nListSize);
//     }
// }

/*
* Neighbourlist update with neighbour cells
*/

void updateNeighbourList()
{
    int count = 0;
    for(int ip=0; ip<np; ip++){
        nPtr[ip] = count;
           
            int iIndex = floor((demPart[ip].pos[0] - xmin) / domainDx);
            int jIndex = floor((demPart[ip].pos[1] - ymin) / domainDy);
            int kIndex = floor((demPart[ip].pos[2] - zmin) / domainDz);
            int cI = iIndex + jIndex * xDiv + kIndex * xDiv * yDiv;
            for (int r = kIndex - 1; r < kIndex + 2; r++)
            {
                for (int q = jIndex - 1; q < jIndex + 2; q++)
                {
                    for (int p = iIndex - 1; p < iIndex + 2; p++)
                    {
                        if(p >= 0 && q >= 0 && r >= 0 && p < xDiv && q < yDiv && r < zDiv)
                        {
                            
                        int neighCellIndex = p + q * xDiv + r * xDiv * yDiv;
                        if (neighCellIndex > xDiv * yDiv * zDiv)
                        {
                            writeLogNum("logfile2.log", "neighCellIndex>xDiv*yDiv*zDiv ", neighCellIndex);
                        }
                        if (iIndex * jIndex * kIndex < 0)
                        {
                            writeLogNum("logfile2.log", "iIndex*jIndex*kIndex < 0 ", iIndex * jIndex * kIndex);
                        }

                        for (int j = 0; j < bdBox[neighCellIndex].noOfParticles; j++)
                        {
                            int jp = bdBox[neighCellIndex].parts[j];
                            vecSub(demPart[ip].pos,demPart[jp].pos,tempDist);
                
                            if(vecMag(tempDist) < 1.55*0.5*(demPart[ip].dia + demPart[jp].dia) && jp > ip)
                            //if (getCenterDist(ip, jp) < rOut && jp > ip)
                            {
                                neighList[count] = jp;
                                count++;
                            }
                        }
                        }
                    }
                }
            }
            
        //delete particle from previous cell
        deleteParticle(ip, demPart[ip].prevCellIndex);
        
        //insert to new cell
        insertToBdBox(ip, cI);
        demPart[ip].displacement[0] = 0.0;
        demPart[ip].displacement[1] = 0.0;
        demPart[ip].displacement[2] = 0.0;
        demPart[ip].bond = 0;
    }
    nListSize = count;
}
