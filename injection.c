
#include "common.h"

void injectParticles(){
    np = 0;
 
    for(int i=0; i<parArraySize; i++)
    {
        //Set particle mass according to DEM density input
        //P_MASS(p) = (4.0/3.0)*PI*pow((0.5*P_DIAM(p)),3.0)*dens;
        smallestParDia = fmin(smallestParDia, demPart[i].dia);
        largestParDia = fmax(largestParDia, demPart[i].dia);
        np++;
    }


    //Setup DEM scaling 
    
    writeLogNum("initial.log","smallestParDia ",smallestParDia);
    writeLogNum("initial.log","largestParDia ",largestParDia);
    writeLogNum("initial.log","minRotEdgeLength ", minRotEdgeLength);
    writeLogNum("initial.log","maxRotEdgeLength ", maxRotEdgeLength);
    writeLogNum("initial.log","minStatEdgeLength ", minEdgeLength);
    writeLogNum("initial.log","maxStatEdgeLength ", maxEdgeLength);
    writeLogLine("initial.log","setReduceUnits()");
    setReduceUnits();
    
    //Assign wall faces to bounding box
    // assignFaceToBDBox(vertX, vertY, vertZ, tempCon1, tempCon2, tempCon3,
    //     noOfVertices,noOfFaces, domainDx, domainDy, domainDz, xDiv, yDiv, 1);
    // newCon1 = allocateIntArray(noOfFacesNew);
    // newCon2 = allocateIntArray(noOfFacesNew);
    // newCon3 = allocateIntArray(noOfFacesNew);
    // writeLogNum("initial.log","noOfFacesNew ",noOfFacesNew);

    // insertToSortedFace(xDiv, yDiv, zDiv,
    //     tempCon1, tempCon2, tempCon3, newCon1, newCon2, newCon3, cellPtr, 1);
    
    writeLogLine("wallmesh.log","assignFaceToBDBox\n");
    writeLogLine("wallmesh.log","assignRotFace\n");

    writeLogLine("wallmesh.log","refineBdBoxFace\n");
    writeLogNum("initial.log","max faces in a BD box (stator) ",maxFaces);
    writeLogNum("initial.log","max faces in a BD box (rotor) ",maxRotFaces);


    //Insert particles to cell 
 
    for(int i=0; i<parArraySize; i++)
    {
                //Insert to cell
                
                demPart[i].mass = (4.0/3.0)*PI*pow((0.5*1e-3*demPart[i].dia),3.0)*dens*massFactor;
                demPart[i].pos[0] =  demPart[i].pos[0]*lengthFactor*1e-3;
                demPart[i].pos[1] =  demPart[i].pos[1]*lengthFactor*1e-3;
                demPart[i].pos[2] =  demPart[i].pos[2]*lengthFactor*1e-3;
                demPart[i].vel[0] = demPart[i].vel[0]*velocityFactor;
                demPart[i].vel[1] = demPart[i].vel[1]*velocityFactor;
                demPart[i].vel[2] = demPart[i].vel[2]*velocityFactor;
                demPart[i].dia = demPart[i].dia*lengthFactor*1e-3;
                
                
                demPart[i].inert = 2.0*demPart[i].mass*pow(0.5*demPart[i].dia,2)/5.0; 
                //initial forces
                demPart[i].force[0] = 0.0;//demPart[i].mass;; //gravitational force
                demPart[i].force[1] = 0.0;//-demPart[i].mass;//gravitational force
                demPart[i].force[2] = -demPart[i].mass; 
                demPart[i].momentum[0] = 0.0;
                demPart[i].momentum[1] = 0.0;
                demPart[i].momentum[2] = 0.0;  
    }
    
    int maxBDSize = 0;
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        if(bdBox[i].noOfParticles > maxBDSize){
        maxBDSize = bdBox[i].noOfParticles;
        }
    }
    writeLogNum("initial.log","maxBDSize ",maxBDSize);

    int maxNBSize = 0;
    for(int i=0; i<np; i++){
        if(maxNBSize, demPart[i].noOfNeigh > maxNBSize){
        maxNBSize = demPart[i].noOfNeigh;
        }
    }

    //initialized = 1;
    writeLogNum("initial.log","maxNBSize ",maxNBSize);

    writeLogNum("initial.log","XMIN ",xmin/lengthFactor);
    writeLogNum("initial.log","YMIN ",ymin/lengthFactor);
    writeLogNum("initial.log","ZMIN ",zmin/lengthFactor);
    writeLogLine("initial.log","----------------\n");
    writeLogNum("initial.log","XMAX ",xmax/lengthFactor);
    writeLogNum("initial.log","YMAX ",ymax/lengthFactor);
    writeLogNum("initial.log","ZMAX ",zmax/lengthFactor);
    writeLogLine("initial.log","----------------\n");
    writeLogNum("initial.log","X LENGTH ",(xmax-xmin)/lengthFactor);
    writeLogNum("initial.log","Y LENGTH ",(ymax-ymin)/lengthFactor);
    writeLogNum("initial.log","Z LENGTH ",(zmax-zmin)/lengthFactor);
    writeLogLine("initial.log","----------------\n");
    writeLogNum("initial.log","DOMAIN DX ",domainDx/lengthFactor);
    writeLogNum("initial.log","DOMAIN DY ",domainDy/lengthFactor);
    writeLogNum("initial.log","DOMAIN DZ ",domainDz/lengthFactor);
    writeLogLine("initial.log","----------------\n");
    writeLogNum("initial.log","DOMAIN XDIV ",xDiv);
    writeLogNum("initial.log","DOMAIN YDIV ",yDiv);
    writeLogNum("initial.log","DOMAIN ZDIV ",zDiv);
    //writeLogNum("initial.log","NUM OF WALL FACES ",noOfFaces)
    //writeLog3Num("initial.log","wall grid divisions ",wBXDiv,wBYDiv,wBZDiv);
    //writeLog3Num("initial.log","wall grid size ",wBoxDx/lengthFactor,wBoxDy/lengthFactor,wBoxDz/lengthFactor);
    writeLogNum("initial.log","MIN EDGE LENGTH ",minEdgeLength);
    writeLogNum("initial.log","MAX EDGE LENGTH ",maxEdgeLength);
    writeLogNum("initial.log","max faces in a BD box ",maxFaces);
    writeLogNum("initial.log","no of rotor faces ",noOfRotorFaces);
    writeLogNum("initial.log","no of stator faces ",noOfStatorFaces);
    writeLogLine("initial.log","----------------\n");
}