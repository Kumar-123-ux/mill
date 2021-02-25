#include "common.h"


/*Read particle information and domain information
*/

void wallInit()
{
  initializePtr();
  
  xmin = 1.e6;
  ymin = 1.e6;
  zmin = 1.e6;
  xmax = -1.e6;
  ymax = -1.e6;
  zmax = -1.e6;  
  demRotAng = 0.;
  direction = 1.0;
  initialized = 0;
  noOfCFDCells = 0;
  minEdgeLength = 1.0e3;
  maxEdgeLength = 0.;
  minRotEdgeLength = 1.0e3;
  maxRotEdgeLength = 0.;
  minFaceArea = 1.e3; 
  maxFaceArea = 0.;
  noOfRotorFaces = 0;
  torqueRadius = 0.0;
  totalTorque = 0.0;
  rotAng = 0.;
  axis1 = 0; // zaxis
  axis2 = 2; // zaxis

  noOfNodes = 0;
  noOfFaces = 0;
  sim_time = 0.0;
  rotorAngPosition = 0.0;

  int noOfThreads = 0;
  int n;

  getBoxSize("infile");
  writeLogNum("initial.log","BOXSIZE ",boxsize);
 
//cfdcell = allocateCFDCell(noOfCFDCells);
//   for(int i=0; i<noOfCFDCells; i++){
//         cfdcell[i].porosity = 0.99; //initial cfdcell porosity is set to 1
//         cfdcell[i].solidVol = 0.0;
//         cfdcell[i].dragFX = 0.0;
//         cfdcell[i].dragFY = 0.0;
//         cfdcell[i].dragFZ = 0.0;
//         cfdcell[i].noOfParts = 0;
//   }
}


void demInit(){
    deallocate();
    //impactExist = false;
    smallestParDia = 1.e3;
    largestParDia = -1.0;   
    part_escaped = 0;
    wall_escaped = 0;
    impAng = 0.0;
    impVel = 0.0; 
    totalImpacts = 0;
    ppTotalImpacts = 0;
    rebldNList = true;

    statZoneID = 3;
    maxFaces = 0;
//    maxRotFaces = 0;
    noOfFacesNew = 0;
//    noOfFacesNewR = 0;
//    noOfRotorFaces = 1;
    noOfStatorFaces = 1;
    maxCnt = 0;
    maxCharge = 0.0;
    maxES = 0.0;
    maxVF = 0.0;
    saveDEM = 0;
    updateDPM = 0;
    //dump_count = 0;
    cfdcycles = 0;
    demcycles = 0;
    //time_count = 0;// time counter used for saving output file
    particle_counter = 0; //counter for keeping the track of number of particles
    demTime = 0.0;
    rotAng = 0.0;
    nListSize = 0;

    readDomain("infile");

    if(!readNp("sample")){
        printf("Particle number cannot read\n");
    }

    //Read material data
    readInput("infile", &dens, &ymod, &pois, &sfc, &rec, &dmpn, &rf, &cyldia, &timeStep, &maxTime,
            &noOfWalls, &updateDPM, &maxVel);
    
    //Read particle-wall contact surfaces
    readWalls("infile", walls);



   //Get wall mesh
    // getWallMeshSize("stator-wall.mesh", &noOfFaces, &noOfVertices);
    // writeLogNum("initial.log","noOfFaces stator ",noOfFaces);
    // writeLogNum("initial.log","noOfVertices stator ",noOfVertices);
    // tempCon1 = allocateIntArray(noOfFaces);
    // tempCon2 = allocateIntArray(noOfFaces);
    // tempCon3 = allocateIntArray(noOfFaces);
    
    // vertX = allocateDoubleArray(noOfVertices);
    // vertY = allocateDoubleArray(noOfVertices);
    // vertZ = allocateDoubleArray(noOfVertices);

    // readConnectivity("stator-wall.mesh", tempCon1, tempCon2, tempCon3);
    // readVertices("stator-wall.mesh", vertX, vertY, vertZ);
    
    // getMinMaxFaceArea(noOfFaces, vertX, vertY, vertZ,
    //   tempCon1, tempCon2, tempCon3, 1);

    //Set min, max values of simulation domain
    xmin = -0.5*cylR;
    xmax = 0.5*cylR;
    ymin = 0.0;
    ymax = cylL;
    zmin = -0.5*cylR;
    zmax = 0.5*cylR;


    //Initial cell divisions for particles
    writeLog3Num("initial.log", "Initial min ",xmin,ymin,zmin);
    writeLog3Num("initial.log", "Initial max ",xmax,ymax,zmax);

    xDiv = floor((xmax-xmin)/boxsize);
    yDiv = floor((ymax-ymin)/boxsize);
    zDiv = floor((zmax-zmin)/boxsize);

    writeLogNum("initial.log","Initial xDiv ",xDiv);
    writeLogNum("initial.log","Initial yDiv ",yDiv);
    writeLogNum("initial.log","Initial zDiv ",zDiv);
    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    //Adjust boundary
    xmin += -domainDx*3.;
    ymin += -domainDy*3.;
    zmin += -domainDz*3.;

    xmax += domainDx*3.;
    ymax += domainDy*3.;
    zmax += domainDz*3.;

    //New cell divisions
    //Initial cell divisions for particles
    xDiv = floor((xmax-xmin)/boxsize);
    yDiv = floor((ymax-ymin)/boxsize);
    zDiv = floor((zmax-zmin)/boxsize);

    writeLogNum("initial.log","FINAL xDiv ",xDiv);
    writeLogNum("initial.log","FINAL yDiv ",yDiv);
    writeLogNum("initial.log","FINAL zDiv ",zDiv);
    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;
    writeLogNum("initial.log","FINAL domainDx ",domainDx);
    writeLogNum("initial.log","FINAL domainDy ",domainDy);
    writeLogNum("initial.log","FINAL domainDz ",domainDz);

    writeLogNum("initial.log","FINAL XMIN ",xmin);
    writeLogNum("initial.log","FINAL YMIN ",ymin);
    writeLogNum("initial.log","FINAL ZMIN ",zmin);
    
    writeLogNum("initial.log","FINAL XMAX ",xmax);
    writeLogNum("initial.log","FINAL YMAX ",ymax);
    writeLogNum("initial.log","FINAL ZMAX ",zmax);

    writeLogNum("initial.log", "No of cells ", noOfCFDCells);
    writeLogNum("initial.log", "No of faces ", noOfFaces);
    writeLogNum("initial.log", "No of nodes ", noOfNodes);
    
    writeLogNum("initial.log", "rotAng (degress) ",rotAng*180./PI);

    writeLogNum("initial.log","PP haConst (*1e20) ",haPP*1e20);
    writeLogNum("initial.log","PP sliding friction ",sfc);
    writeLogNum("initial.log","PP damping ",dmpn);
    writeLogNum("initial.log","PP rolling friction ",rf);
    
    writeLogNum("initial.log","PW haConst (*1e20) ",haPW*1e20);
    writeLogNum("initial.log","PW sliding friction ",pwsf);
    writeLogNum("initial.log","PW damping ",pwdmpn);
    
    writeLogNum("initial.log","lamda1 ",lamda1*1e9);
    writeLogNum("initial.log","lamda2 ",lamda2*1e9);
    writeLogNum("initial.log","rms1 ",rms1*1e9);
    writeLogNum("initial.log","rms2 ",rms2*1e9);
    writeLogNum("initial.log","Update DPM ",updateDPM);
    writeLogNum("initial.log","Particle Array Size ",parArraySize);
    writeLogNum("initial.log","density ",dens);
    writeLogNum("initial.log","Youngs Modulus ",ymod);
    writeLogNum("initial.log","largestParDia ",largestParDia);
    writeLogNum("initial.log","smallestParDia ",smallestParDia);
    writeLogNum("initial.log","V1 ",V1);
    writeLogNum("initial.log","V2 ",V2);
    writeLogNum("initial.log","imageConst ",imageConst);
    writeLogNum("initial.log","Zs (x1e9) ",Zs*1.e9);
    writeLogNum("initial.log","alpha ",alpha);
    writeLogNum("initial.log","permitivity (x1e-12)",permitivity*1e12);
    writeLogNum("initial.log","K0 ",imageConst*Zs/permitivity/(4.0*PI*pow(largestParDia*0.5,2)));
    writeLogNum("initial.log","esfTrue ",esfTrue);
    writeLogNum("initial.log","capfTrue ",capfTrue);
    double qMax = (V1-V2)/(imageConst*Zs/permitivity/(4.0*PI*pow(largestParDia*0.5,2)));
    double maxESForce = qMax*qMax/(4.*PI*permitivity*1.0e-6*1.0e-6);
    writeLogNum("initial.log","Qmax (x1e12)",qMax*1e12);
    writeLogNum("initial.log","maxESForce (x1e6) ",maxESForce*1e6);
    writeLogNum("initial.log","Timestep (x1e6)",timeStep*1.e6);
    writeLogNum("initial.log","noOfCFDCells ",noOfCFDCells);
    writeLogNum("initial.log","maxBound ",maxBound);
    writeLogNum("initial.log","rotor zone id ",rotZoneID);
    writeLogNum("initial.log","rotor angular velocity ",angVel);
    writeLogNum("initial.log","Min separation  ",s_min);
    writeLogNum("initial.log","Liquid volume (x1e9) ",liq_vol*1e9);
    writeLogNum("initial.log","Surface tension  ",surf_tens);
    writeLogNum("initial.log","Contact angle  ",cont_ang);

     //Allocate memory
    allocate();
    writeLogLine("logfile.log"," allocate\n ");
    writeLogNum("initial.log","MIN EDGE LENGTH (mm)",minEdgeLength*1.e3);
    writeLogNum("initial.log","MAX EDGE LENGTH (mm)",maxEdgeLength*1.e3);
    writeLogNum("initial.log","MIN FACE AREA (x1.e12) ",minFaceArea*1.e12);
    writeLogNum("initial.log","MAX FACE AREA (x1.e12) ",maxFaceArea*1.e12);
    writeLogNum("initial.log","Drum Radius ",0.5*cylR);
    writeLogNum("initial.log","Drum Length ",cylL);
    writeLogNum("initial.log","Lifter Width ",lifterWidth);
    writeLogNum("initial.log","Lifter Height ",lifterHeight);
    
    // Read particle position
    if(!readParticles("sample")){
        // for (int i=0; i<parArraySize; i++)
        // {
        //     demPart[i].pos[0] = (demPart[i].pos[0])*1e-3;
        //     demPart[i].pos[1] = (demPart[i].pos[1])*1e-3;
        //     demPart[i].pos[2] = (demPart[i].pos[2])*1e-3;
        //     demPart[i].dia = (demPart[i].dia)*1e-3;
        // }
        printf("Particles cannot read\n");
    }

    writeLogNum("initial.log","SIMTIME ",sim_time);
}

/*
Add wall face to BD box
fc - face index
type - face type (rotor=1, stator=2)
cI - BD box index
*/

void addFaceToBDBox(int fc, int cI, int type){
    if(type == 1){
        bdBox[cI].face[bdBox[cI].totalFaces] = fc;
        bdBox[cI].totalFaces++;
        noOfFacesNew++; //total stator faces
    }
    // else{
    //     bdBox[cI].faceR[bdBox[cI].totalFacesR] = fc;
    //     bdBox[cI].totalFacesR++;
    //     noOfFacesNewR++; //total rotor faces
    // }
}

/* Deallocate memory*/
void deallocate(){
    if(sizeof(walls) != 0){
        free(impactPart);
        free(impactPartPP1);
        free(impactPartPP2);
        free(impactVel);
        free(impactAng);
        free(parDia);
        free(impactPosX);
        free(impactPosY);
        free(impactPosZ);
        free(impactSurface);
        
        free(ppImpactVel);
        free(ppParDia1);
        free(ppParDia2);
        free(ppImpactAng);
        free(ppImpactPosX);
        free(ppImpactPosY);
        free(ppImpactPosZ);

        free(walls);
        free(uVec);
        free(ipRVec);
        free(jpRVec);
        free(torque);
        free(bdBox);
        free(ijVec);
        free(rotVel);
        free(ipCntPntVel);
        free(jpCntPntVel);
        free(cntPntVel);
        free(demPart);
        free(holePart);
        free(inipos);
        //free(wFace);
        //free(wFaceRot);       
        free(vertX);
        free(vertY);
        free(vertZ);
        free(cellPtr);
        // free(vertXR);
        // free(vertYR);
        // free(vertZR);
        // free(cellPtrR);
        if(tempCon1 != NULL){
            free(tempCon1);
            free(tempCon2);
            free(tempCon3);
        }
        // if(tempConR1 != NULL){
        //     free(tempConR1);
        //     free(tempConR2);
        //     free(tempConR3);
        // }
    } 
}

/* Allocate arrays */
void allocate(){
    //If allocated then dealllocate
    //deallocate();
    
    neighList = allocateIntArray(NEIGH_LIST_ARRAY_SZ);
    for(int i=0; i<NEIGH_LIST_ARRAY_SZ; i++){
        neighList[i] = 0;
    }
    nPtr = allocateIntArray(parArraySize);
    for(int i=0; i<parArraySize; i++){
        nPtr[i] = 0;
    }

    // holeCenter1 = allocateDoubleArray(DIM);
    // holeCenter2 = allocateDoubleArray(DIM);
    // for (int i=0; i<DIM; i++){
    //   holeCenter1[i] = 0.0;
    //   holeCenter2[i] = 0.0;
    // }
    printf("%d\n",parArraySize);
    int expand = 10000;
    datacount = 0;
    impactPartPP1 = allocateIntArray(parArraySize + expand);
    impactPartPP2 = allocateIntArray(parArraySize + expand);
    impactPart = allocateIntArray(parArraySize + expand);
    impactVel = allocateDoubleArray(parArraySize + expand);
    impactAng = allocateDoubleArray(parArraySize + expand);
    parDia = allocateDoubleArray(parArraySize + expand);
    impactSurface = allocateIntArray(parArraySize + expand);
    impactPosX = allocateDoubleArray(parArraySize + expand);
    impactPosY = allocateDoubleArray(parArraySize + expand);
    impactPosZ = allocateDoubleArray(parArraySize + expand);
      
    ppImpactVel = allocateDoubleArray(parArraySize + expand);
    ppImpactAng = allocateDoubleArray(parArraySize + expand);
    ppParDia1 = allocateDoubleArray(parArraySize + expand);
    ppParDia2 = allocateDoubleArray(parArraySize + expand);
    ppImpactPosX = allocateDoubleArray(parArraySize + expand);
    ppImpactPosY = allocateDoubleArray(parArraySize + expand);
    ppImpactPosZ = allocateDoubleArray(parArraySize + expand);

    uVec = allocateDoubleArray(DIM);
    ipRVec = allocateDoubleArray(DIM);
    jpRVec = allocateDoubleArray(DIM);
    ijVec = allocateDoubleArray(DIM);
    rotVel = allocateDoubleArray(DIM);
    torque = allocateDoubleArray(DIM);
    ipCntPntVel = allocateDoubleArray(DIM);
    jpCntPntVel = allocateDoubleArray(DIM);
    cntPntVel = allocateDoubleArray(DIM);
  
    bdBox = allocateBdBoxArray(xDiv*yDiv*zDiv); //bounding box array
    //cellPtr = allocateIntArray(xDiv*yDiv*zDiv); //cell array index poinitng to sorted face array
    //cellPtrR = allocateIntArray(xDiv*yDiv*zDiv); //cell array index poinitng to sorted rotor face array

    writeLogLine("initial.log"," allocate stage 1\n ");

    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        bdBox[i].noOfParticles = 0;
        bdBox[i].totalFaces = 0;
        //bdBox[i].totalFacesR = 0;
        bdBox[i].face = allocateIntArray(MAX_FACES);
        bdBox[i].parts = allocateIntArray(NO_OF_PARTICLES_IN_BDCELL);
        //bdBox[i].faceR = allocateIntArray(MAX_FACES);
    }
   
    //tempContactList = allocateIntArray(parArraySize);
    demPart = allocatePar(parArraySize);
    holePart = allocateIntArray(parArraySize);
    inipos = allocateIniData(parArraySize);

    for(int i=0; i<parArraySize; i++){
        demPart[i].injtime = 0.0;
        demPart[i].rotMomc = 0.0;
        demPart[i].trapped = false;
        holePart[i] = 0;
        demPart[i].ppNo = 1;
        demPart[i].pwNo = 1;
        demPart[i].ppHoleTngF = 0.0;
        demPart[i].ppHoleNrmF = 0.0;
        demPart[i].pwHoleTngF = 0.0;
        demPart[i].pwHoleNrmF = 0.0;
        demPart[i].holeDragF = 0.0;
        demPart[i].holeVelX = 0.0;
        demPart[i].holeVelY = 0.0;
        demPart[i].holeVelZ = 0.0;
        demPart[i].wallCnt = 0;
        demPart[i].elecForce = 0.0;
        demPart[i].capForce = 0.0;
        demPart[i].vanForce = 0.0;
        demPart[i].preRelVelWall = 0.0;
        demPart[i].preRelAngWall = 0.0;
        demPart[i].impactExist = false;
        demPart[i].minWallDist = 1.0e3;
        demPart[i].bond = 0;
        demPart[i].nrmDisp = 0.0;
        // demPart[i].impactVel = 0.0;
        // demPart[i].impactAng = 0.0;
        //demPart[i].noOfNeigh = 0;
        demPart[i].cordNo = 0;
        demPart[i].pos = allocateDoubleArray(DIM);
        demPart[i].angVel = allocateDoubleArray(DIM);
        demPart[i].vel = allocateDoubleArray(DIM);
        demPart[i].fVel = allocateDoubleArray(DIM);
        //demPart[i].hisDisp = allocateDoubleArray(DIM);

        // demPart[i].hisDisp[0]  = 0.0;
        // demPart[i].hisDisp[1]  = 0.0;
        // demPart[i].hisDisp[2]  = 0.0;   

        demPart[i].angVel[0] = 0.0;
        demPart[i].angVel[1] = 0.0;
        demPart[i].angVel[2] = 0.0;

        demPart[i].force = allocateDoubleArray(DIM);
        // demPart[i].dragForce = allocateDoubleArray(DIM);
        demPart[i].momentum = allocateDoubleArray(DIM);
        demPart[i].momentum[0] = 0.0;
        demPart[i].momentum[1] = 0.0;
        demPart[i].momentum[2] = 0.0;
        // demPart[i].dragForce[0] = 0.0;
        // demPart[i].dragForce[1] = 0.0;
        // demPart[i].dragForce[2] = 0.0;
        for(int j=0; j<DIM; j++){
            demPart[i].displacement[j] = 0.0;
        }
        demPart[i].insertable = 1;
        demPart[i].active = false;
        demPart[i].escaped = false;
        demPart[i].inserted = false;
        
        demPart[i].noOfCnt = 0;
        demPart[i].noOfWallCnt = 0;
        demPart[i].incontact = 0;
        demPart[i].eCharge = 0.0;//1.03e-4*pow(0.5*largestParDia,1.7);
        demPart[i].noOfPartColl = 0;
        demPart[i].noOfWallColl = 0;
        demPart[i].totalParColl = 0;

        for(int j=0; j<CNTSIZE; j++){
            demPart[i].contList[j] = 0;
            demPart[i].impVelPre[j] = 0.0;
            demPart[i].impAngPre[j] = 0.0;
            demPart[i].impactExistPP[j] = false;
            //demPart[i].bondList[j] = 0;
            
        }
        // for(int j=0; j<WALLCNTSIZE; j++){
        //     demPart[i].wallCntList[j] = 0;
        // }
	    //demPart[i].maxPartCollE = 0.0;
        //demPart[i].maxWallCollE = 0.0;
        //demPart[i].preFluidCell = -1;
        //demPart[i].wallType = 0;
    }

    walls = allocateIntArray(noOfWalls);
    writeLogLine("initial.log"," allocate completed\n ");
}

/*
Assign scale factors
*/
void setReduceUnits()
{
    //Scale factors for reduced units
	refLength = largestParDia;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravity/refLength);
	densityFactor = 6.0/(PI*refDensity);
	forceFactor = 6.0/(gravity*PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravity*PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravity*PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((double)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = 1.2*largestParDia*lengthFactor;
    contactGap = 0.01; //contact gap limit as a percentage of particles size
    dsmaxCff = sfc*(2.0-pois)/(2.0*(1.0-pois));
    writeLogNum("initial.log"," mg in reduced units",(4.0/3.0)*PI*pow(largestParDia*0.5,3)*largestParDensity*massFactor);
    
    dti = 0.0;
    dd = 0.0;
    dsmax = 0.0;

    // scale particle properties
    ymod = ymod*pressureFactor;
    cyldia = cyldia*lengthFactor;
	timeStep = timeStep*timeFactor;
	maxTime	= maxTime*timeFactor;
    angVel = (2.0*M_PI*angVel)/(60.0*timeFactor); // rad sec in reduced units
    //haConst = 1.6E-19;
    haPP = haPP*forceFactor*lengthFactor;
    haPW = haPW*forceFactor*lengthFactor;
    dmpn = dmpn*timeFactor;
    //pwdmpn = pwdmpn*forceFactor/(pressureFactor*lengthFactor*velocityFactor);
    pwdmpn = pwdmpn*timeFactor;
    prevCPUTime = 0.;

    elasticMod = ymod/(1.0-pow(pois,2));
   //Find allowed displacment for neighbourlist update
    rIn = smallestParDia*lengthFactor;
    rOut = 1.55*rIn; //By definition
    //allowedDisp = 0.5*(rOut-rIn);
    allowedDisp = 0.5*(rOut-rIn);
    s_rupture = (1+0.5*cont_ang)*pow(liq_vol,1/3);

    for (int i=0; i<np; i++){
        demPart[i].haPp = haPP;
        //demPart[i].haPw = haPW;
        //demPart[i].dt = timeStep;
        for(int j=0; j<DIM; j++){
            demPart[i].displacement[j] = 2.0*allowedDisp;
        }
    }

    //Multiply node coordinates by lengthFactor
    // for(int i=0; i<noOfNodes; i++)
    // {
    //     ndArray[i].coord[0] = ndArray[i].coord[0]*lengthFactor;
    //     ndArray[i].coord[1] = ndArray[i].coord[1]*lengthFactor;
    //     ndArray[i].coord[2] = ndArray[i].coord[2]*lengthFactor;        
    // }

    //Multiply face nodes by lengthFactor
    //initializeFaces(vertX, vertY, vertZ, noOfVertices);
    //initializeFaces(vertXR, vertYR, vertZR, noOfVerticesR); 
    //Adjust boundary limits
    xmin = xmin*lengthFactor;
    ymin = ymin*lengthFactor;
    zmin = zmin*lengthFactor;
    xmax = xmax*lengthFactor;
    ymax = ymax*lengthFactor;
    zmax = zmax*lengthFactor;

    domainDx = domainDx*lengthFactor;
    domainDy = domainDy*lengthFactor;
    domainDz = domainDz*lengthFactor;

    cellRadius = 0.5*sqrt(domainDx*domainDx+domainDy*domainDy);

    //iter = CURRENT_TIMESTEP/(timeStep/timeFactor);
    //set current time
    //demTime = CURRENT_TIMESTEP*cfdcycles*timeFactor;
    //current rot angle
    //rotAng = -CURRENT_TIMESTEP*cfdcycles*angVel;

    lamda1 = lamda1*lengthFactor;
    lamda2 = lamda2*lengthFactor;
    rms1 = rms1*lengthFactor;
    rms2 = rms2*lengthFactor;
    double k0 = imageConst*Zs/permitivity/(4.0*PI*pow(0.5*largestParDia,2));
    cylR = 0.5*lengthFactor*cylR; //radius
    cylL = lengthFactor*cylL;
    lifterWidth = lengthFactor*lifterWidth;
    lifterHeight = lengthFactor*lifterHeight;
    // holeCenter1[0] = 0.5*cylL+cylR;
    // holeCenter1[1] = 0.0;
    // holeCenter1[2] = 0.0;
    // holeCenter2[0] = -(0.5*cylL+cylR);
    // holeCenter2[1] = 0.0;
    // holeCenter2[2] = 0.0;
    //rf = rf*lengthFactor;
    writeLogNum("initial.log","Start time ",demTime/timeFactor);
    writeLogNum("initial.log","rolling friction ",rf);
    writeLogNum("initial.log","max charge (nC)",1e9*(V2-V1)/k0);
    //writeLogNum("initial.log","DEM iteration CURRENT_TIMESTEP/(timeStep/timeFactor) ",iter);
}

/* Multiply wall face coordinate by lenghFactors*/
void initializeFaces(double *vX, double *vY, double *vZ, int size){
    for(int i=0; i<size; i++){
        vX[i] = vX[i]*lengthFactor;
        vY[i] = vY[i]*lengthFactor;
        vZ[i] = vZ[i]*lengthFactor;
    }
}

/* Resize array */ 
// int *resize_array(int *myint, int sz){
// 	int *tmp = (int*)realoc(myint, sz*sizeof(int));
// 	if(tmp){
// 		myint = tmp;
// 	}
// 	return myint;
// }

/*
Allocate an integer type array
return: int* 
*/
int *allocateIntArray(int size)
{
    int *val = (int*)malloc(size*sizeof(int));
    memset(val,0,size*sizeof(int));
    return val;
}

/*
Allocate a double* type array
return: double* 
*/
double *allocateDoubleArray(int size)
{
    double *val = (double*)malloc(size*sizeof(double));
    memset(val,0.0,size*sizeof(double));
    return val;
}

/*
Allocate a char* type array
return: char* 
*/
char *allocateCharArray(int size)
{
    char *val = (char*)malloc(size*sizeof(char));
    return val;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct BdBox *allocateBdBoxArray(int size)
{
    struct BdBox *bdB = (struct BdBox*)malloc(size*sizeof(struct BdBox));
    return  bdB;
}

/*
Allocate particle array
*/
struct demParticle *allocatePar(int np)
{
    struct demParticle *par = (struct demParticle*)malloc(np*sizeof(struct demParticle));
    return par;
}

/*
Allocate inidata array
*/
struct inidata *allocateIniData(int np)
{
    struct inidata *ps = (struct inidata*)malloc(np*sizeof(struct inidata));
    return ps;
}

/*
Allocate wall face array
*/
struct wallFace *allocateFace(int nc)
{
    struct wallFace *fc = (struct wallFace*)malloc(nc*sizeof(struct wallFace));
    return fc;
}

/*
Allocate wall node array
*/
struct wallNode *allocateNode(int nn)
{
    struct wallNode *nd = (struct wallNode*)malloc(nn*sizeof(struct wallNode));
    return nd;
}

/*
Allocate cfd eell array
*/
struct cfdCell *allocateCFDCell(int nc)
{
    struct cfdCell *cfdcell = (struct cfdCell*)malloc(nc*sizeof(struct cfdCell));
    return cfdcell;
}

/*
* Don't forget to initialize pointers
*/
void initializePtr(){
    tempCon1 = NULL;
    tempCon2 = NULL;
    tempCon3 = NULL;
    newCon1 = NULL;
    newCon2 = NULL;
    newCon3 = NULL;
    vertX = NULL;
    vertY = NULL;
    vertZ = NULL;
    // tempConR1 = NULL;
    // tempConR2 = NULL;
    // tempConR3 = NULL;
    // newConR1 = NULL;
    // newConR2 = NULL;
    // newConR3 = NULL;
    // vertXR = NULL;
    // vertYR = NULL;
    // vertZR = NULL;
    cellPtr = NULL;
    //cellPtrR = NULL;

    nPtr = NULL;
    neighList = NULL;
    
    walls = NULL;
    uVec = NULL;
    ipRVec = NULL;
    jpRVec = NULL;
    torque = NULL;
    ijVec = NULL;
    rotVel = NULL;
    ipCntPntVel = NULL;
    jpCntPntVel = NULL;
    cntPntVel = NULL;
}