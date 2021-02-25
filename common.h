
#ifndef COMMON_H
#define COMMON_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

/*--------------Common definitions---------------*/
//#define NUM_MAT 2    // two types material
#define PI 3.1415926f
#define gravity 9.8f //ms^-2
#define conversion 1.0e-3f //convert length values to meters
#define overlapLimit 0.0005f //ovelap cuttoff value

/*------------ Particle information -------------*/
#define XMAX_BOUND 0.015 //Maxmimum x bound for active particles
//#define multif3 20//multification factor used in bounding box division
#define DIM 3 // 3D problem
//#define NBSIZE 1 //size of neighbourlist
#define CNTSIZE 40
#define MAX_FACES 10
//#define WALLCNTSIZE 3 //Maximum number of walls, a particle can have contact with
//#define NO_OF_FACES 4//number of faces contacting with particles
#define NO_OF_PARTICLES_IN_BDCELL 40
#define NEIGH_LIST_ARRAY_SZ 1000000
#define MAXIMPACTS 1000

// err detection utility
#define FPRINTF(a) fprintf a

/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/

//FILE *LogFile;
bool rebldNList;
double demRotAng;
short int direction;
int ftime;
double minEdgeLength,maxEdgeLength, maxBound;
double minRotEdgeLength, maxRotEdgeLength;
double minFaceArea, maxFaceArea;
double prevCPUTime;
int cfdcycles, demcycles, iter, noOfIteration;
int updateSourceTerm, rotZoneID, statZoneID;
int particle_counter;
int noOfWalls, noOfNodes, noOfFaces, noOfCFDCells, noOfRotorFaces, noOfStatorFaces;
short int maxFaces, maxRotFaces, nFace, faceType; //temp value

int *walls;
int np, parArraySize; //number of particles, particle array size
double sim_time, rotorAngPosition, torqueRadius, totalTorque; 
unsigned int initialized,axis1, axis2;
 
double inletVel, angVel, rotAng; //angVel is the rotor speed in rad/sec
double ductxmin, ductxmax, ductxedge1, ductxedge2, ductymin, ductymax, ductzmin, ductzmax, ductzedge;
double cutGap; //particle parameters
double contactGap; //limit for contact gap
double cellRadius; //Radius defined by bounding box cell given by 0.5*sqrt(dx^2+dy*2)
double refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor, 
	smallestParDia, largestParDia,largestParDensity;
double vGapMn; //minimum gap for vanderwal force 
double maxVel, impVel, impAng, nrmImpVel; //maximum flow velocity at the inlet 
double lamda1, lamda2, rms1, rms2;
double s_min,liq_vol, surf_tens, cont_ang; 
int xDiv, yDiv, zDiv; //Number of divisions in orthogonal domain cells
//int wBXDiv, wBYDiv, wBZDiv; //Number of divisions for wall mesh
double boxsize, domainDx, domainDy, domainDz; //Domain cell size
//double wBoxDx, wBoxDy, wBoxDz; //wall grid cell size 
double xmax, ymax, zmax; //Max values of domain boundary 
double xmin, ymin, zmin; //Min values of domain boundary 

double Zs, V1, V2, permitivity, imageConst, alpha, ks, esfTrue, capfTrue; //electrostatic constants
double rIn, rOut, allowedDisp; //if particle displacement > allowedDisp -> update neighbourlist

double *uVec, *ipRVec, *jpRVec, *ijVec, *rotVel, *torque;
double *ipCntPntVel, *jpCntPntVel, *cntPntVel;
double dens, ymod, pois, sfc, rf, rec, dmpn, elasticMod, haPP, haPW; //particle material property
double pwsf, pwdmpn, s_rupture;
double dispx, dispy, dispz, dsmaxCff, dti, dd, dsmax, fdt; //used in contact force calculation
double fdtx, fdty, fdtz;
double cyldia; //cylinder diameter
double timeStep, demTime, maxTime;
struct BdBox *bdBox, *tbdBox;
int *cellPtr, *cellPtrR;

//struct wallBox *wBox;
struct demParticle *demPart;
struct inidata *inipos;
struct cfdCell *cfdcell;
int *holePart;
//double *impactVel;
int  ppTotalImpacts, datacount;

int updateDPM;
int saveDEM; //counter used for saving DEM particle for TECPLOT
int maxCnt;
double maxCharge, maxES, maxVF;
unsigned int maxChargePart;

int *tempCon1, *tempCon2, *tempCon3;//connectivity array for wall faces
int *newCon1, *newCon2, *newCon3;//connectivity array for wall faces
double *vertX, *vertY, *vertZ;// wall vertices
int noOfVertices, noOfFacesNew;

// int *tempConR1, *tempConR2, *tempConR3;//temporary connectivity array for rotor wall faces
// int *newConR1, *newConR2, *newConR3;//connectivity array for rotor wall faces
// double *vertXR, *vertYR, *vertZR;// rotor wall vertices
// int noOfVerticesR, noOfFacesR, noOfFacesNewR; 

int *nPtr;
int *neighList;
int nListSize;

double cylR, cylL, lifterWidth, lifterHeight;
//double *holeCenter1, *holeCenter2;
double rVec[DIM], cntPnt[DIM], tempUV[DIM], edgeCent[DIM], edgeCnt[DIM], capCent[DIM], cntPntEdge[DIM];
double tngUVec[DIM], fdtVec[DIM], totalForce[DIM], rotMom[DIM], momentum[DIM];
double rotOmega[DIM];
double relativeVel[DIM], tempVel[DIM];

double *impactVel, *parDia, *impactAng, *impactPosX, *impactPosY, *impactPosZ;
double tempDist[DIM], temVec[DIM];
int *impactSurface, *impactPart, *impactPartPP1, *impactPartPP2;
int totalImpacts;


double *ppImpactVel, *ppParDia1, *ppParDia2, *ppImpactAng, *ppImpactPosX, 
		*ppImpactPosY, *ppImpactPosZ;

int part_escaped, wall_escaped;
struct wallNode *ndArray;
struct wallFace *fcArray;
//int *tempContactList;

struct wallNode{
	int nd;
	double coord[DIM];
};

struct wallFace{
	bool rot;
	int wNode[DIM];
};

//cfdcell information
struct cfdCell{
	double porosity;
	double solidVol;
	double dragFX;
	double dragFY;
	double dragFZ;
	int noOfParts;
};


//Bounding box which holds CFD cells and DEM particles
struct BdBox{
	short int maxSize;
	int totalFaces, totalFacesR; //number of wall faces
	int noOfParticles; //number of DEM particles
	int *parts;
	int *face;
	// int *faceR;
};


//Particle
struct demParticle{
	bool trapped;
	double minWallDist, rotMomc, injtime;
	int faceId, faceType, bond;
	double preRelVelWall, preRelAngWall;
	//double impactVel, impactAng; //impactAng is the angle between imapct velocity and unit vector
	double displacement[DIM]; //if displacement > rMax update neighbourlist
	double dia, inert, mass, nrmDisp;
	double *pos, *angVel, *vel, *force, *momentum, *fVel;
	double vanForce, capForce, elecForce;
	int contList[CNTSIZE];
	double impVelPre[CNTSIZE];
	double impAngPre[CNTSIZE];
	bool impactExistPP[CNTSIZE];
	bool impactExist; //wall
	//int wallCntList[WALLCNTSIZE];
	//int bondList[CNTSIZE];
	int noOfCnt,noOfPartColl, noOfWallColl, noOfWallCnt, totalParColl, wallCnt;
	//double maxPartCollE, maxWallCollE;
	int incontact;
	int noOfNeigh, cordNo;
	int prevCellIndex;	
	short int insertable;
	bool active, escaped, inserted; 
	double haPp;//Hamarker constant
	//double *dragForce;
	double dragFX, dragFY, dragFZ;
	double eCharge;
	double ppHoleTngF, ppHoleNrmF, pwHoleTngF, pwHoleNrmF, holeDragF, holeVelX, holeVelY, holeVelZ;
	int pwNo, ppNo, cntType;
};

//Store data for ejecting particles
struct inidata{
	int parno;
	double dt, dia, xpos, ypos, zpos, xvel, yvel, zvel;
};

void wallInit();
void injectParticles();
void updateNeighbourList();
void printCPUTime();
double solidFraction(int ip);
//bool getProjectOnFace(double *p, double *n1, double *n2, double *n3, double *dVec, int ip, int nodeNo);
double triArea(double *v1, double *v2);
void setWallGridSize();
void minMaxEdgeLength(double *n1, double *n2, double *n3, int type);
void getFace(double *n1, double *n2, double *n3, 
   int c1, int c2, int c3, double *vX, double *vY, double *vZ);
void writeDumpFile(char *infile);
void writeInjectionTime(char *inFile);
void readInjectionTime(char *inFile);
void writeInjectionTime(char *inFile);
void writeLogNum(char *infile, char *line, double num);
void writeLogLine(char *infile, char *line);
void writeLog3Num(char *infile, char *line, double v1, double v2, double v3);
void dumpImpact(char *infile);
void dumpPPImpact(char *infile);
void addPPContact(int ip, double impAng, double impVel, int jp);
void addImpact(int ip, double impAng, double impVel, int surface);
double readInputVelocity(char *infile);
void readInput(char *infile, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *rf, double *cyldia,
			 double *dt, double *mT, int *nW, int *updateDPM, double *maxVel);
void readGeom(char *infile, double *ductxmin, double *ductxmax, double *ductxedge1, double *ductxedge2, double *ductymin, 
            double *ductymax, double *ductzmin, double *ductzmax, double *ductzedge);
//void readWallMesh(char *infile, int mode);
//void readFaces(int end, FILE *wFile, struct wallFace *wF);
void readDomain(char *infile);
int readNp(char *infile);
int readParticles(char *infile);
//void test(Tracked_Particle *p, Thread *t);

int *resize_array(int *myint, int sz);
int *allocateIntArray(int size);
double *allocateDoubleArray(int size);

struct tempBdBox *allocateTBdBoxArray(int size);
struct BdBox *allocateBdBoxArray(int size);
struct wallBox *allocateWallBoxArray(int size);
struct wallFace *allocateFace(int nc);
struct wallNode *allocateNode(int nn);
struct demParticle *allocatePar(int np);
struct inidata *allocateIniData(int np);
struct cfdCell *allocateCFDCell(int nc);

double partVol(int p);
void getProjection(double *p, double *n1, double *n2, double *n3, double *dVec);
void getWallMeshSize(char *infile, int *conSize, int *noOfVert);
void readConnectivity(char *infile, int *con1, int *con2, int *con3);
void readVertices(char *infile, double *vertX, double *vertY, double *vertZ);
void getMinMaxFaceArea(int noOfF, double *vX, double *vY, double *vZ,
      int *con1, int *con2, int *con3, int type);
void getFace(double *n1, double *n2, double *n3, 
   int c1, int c2, int c3, double *vX, double *vY, double *vZ);

void insertToBdBox(int p, int cI);
void addToBdBox();
void addFaceToBDBox(int fc, int cI, int type);
void deleteFace(int fc, int type, int cI);

void assignFaceToBDBox(double *vX, double *vY, double *vZ, int *con1, int *con2, int *con3,
						int noOfVert, int noOfF, double domDX, double domDY, double domDZ, int xD, int yD, int type);
void insertToSortedFace(int xDv, int yDv, int zDv,int *con1, int *con2, 
	int *con3, int *nCon1, int *nCon2, int *nCon3, int *cPtr, int type);
void assignToBdBox(struct BdBox *bBox, double *vX, double *vY, double *vZ, 
               int *con1, int *con2, int *con3, 
               int xDv, int yDv, int zDv, double dDx, double dDy, double dDz,
               double edgeLength, int type);

// void readInput(char *infile, int *np, double *dens, double *ymod, 
// 			double *pois, double *sfc, double *rec, double *dmpn, double *rf,
// 			double *cyldia, double *dt, int *nW, int *updateDPM, double *maxVel, int *rotZoneID);
void impactVelocityPP(int ip, int jp, double *uV);
void impactVelocity(int ip, double *uV);
void checkFaceInTriangle(int ip, double *node1, double *node2, double *node3, double rAng, double omega);
void getProjection(double *p, double *n1, double *n2, double *n3, double *dVec);
bool getProjectOnFace(double *p, double *n1, double *n2, double *n3, double *dVec);
void getBoxSize(char *infile);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, struct demParticle *par, int *np);
void readWalls(char *infile, int *walls);
void getNearFace(int ip, int *con1, int *con2, int *con3,
                 double *vX, double *vY, double *vZ, int *cPtr,
                 double dDx, double dDy, double dDz, int xD, int yD, int zD, 
                 int noOfF, double rAng, double omega);
void contactPointInXYCoord(double *cntPnt, int ip, double rAng);

void writeTec();
void demInit();
void buildDEMCFDCellMap();
void copyDEMInfo();

void writeDump();
void demSave();
void writeTorque();
void dumpPPImpact();
void dumpPWImpact();
void writeFluidVelocity();
void cordSave();
void cellSave();
void allocate();
void deallocate();

//void updateForce(Tracked_Particle *p);
void findContactFromMesh(int ip, double *dVec);
//void findContactFromBoundary(Particle *p);
//void calculateDragForce(Particle *p);
void boundaryContactForce(int pI, double *n1, double *n2, double *n3, double *uVec);
void ppVWForce(int ip, int jp, double vGap, double *uV);
void pWallVWForce(int p, double vGap, double *uV);

void assignGravity();
double getOverlap(double *parPos, double dia, double *n1, double *n2, double *n3, double *uVec);

void neighbourContactForce(int pI);
void findDrumContact(int ip, double omega);
void findEndWallContact(int ip, double omega);
void findLifterContact(int ip, double omega);
void surfaceContactForce(int p, double nrmDisp, double *uV);
double getDist(double *v1, double *v2);
void vecAdd(double *v1, double *v2, double *vec);
void crossProd(double *v1, double *v2, double *vec);
void vecSub(double *v1, double *v2, double *vec);
double relVel(int ip, int jp);
void unitVec(double *v, double *vec);
void sclMult(double scl, double *vec);
void sclVecMult(double scl, double *inVec, double *outVec);
double vecMag(double *vec);
void projVec(double *v1, double *n, double *vec, int type);
double dotProduct(double *v1, double *v2);
void getUnitVector(double *v1, double *v2, double *v3, double *uVec);

void initialize(double *sortedList, int *sortedParIndex, int *cellSE, int np,
    double *pos, double *parDia);
void initializeFaces(double *vX, double *vY, double *vZ, int size);
void initializePtr();
////void assignFace(int step, struct wallFace *wF, int zID);
//void refineBdBoxFace(struct wallFace *wF);
void refineBdBox(int xDv, int yDv, int zDv, double dDx, double dDy, double dDz,
   double *vX, double *vY, double *vZ, int *con1, int *con2, int *con3, int type);
int insertable(int ip, int jp);
void addNeighbour(int  ip, int jp);
//void updateNeighbourList(int p);
void deleteNeighbour(int ip, int jp);
void addContact(int ip, int cnt);
void deleteContact(int ip, int cnt);
void deleteAllContacts(int ip, int nF);
void addWallContact(int ip, int cnt);
void deleteWallContact(int ip, int cnt);

int contactExist(int ip, int cnt);
void deleteAllContacts(int ip, int nF);
void update(double *pX, int np);

double getCenterDist(int ip, int jp);
void setReduceUnits();
void updateBBFluidVel();
void insertToBdBox(int p, int cI);
void deleteParticle(int p, int cI);
void forceCalculation(int ip);
void updatePosition(int ip);
void buildNeighList();
void partContactForce(int ip, int jp, double nrmDisp, double *uV);
void surfaceContactForce(int p, double nrmDisp, double *uV);
void dumpPWForce();
void dumpPPForce();

void run();
void test();

#endif 
