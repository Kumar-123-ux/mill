
#include "common.h"


/* Initialize dem information
Executes only once at the begining of FLUEMT*/
int main()
{ 
  //writeLogLine("logfile3.log", "INIT\n");
  wallInit();
  demInit();
  injectParticles();
  // //test();
  run();
  printf("WORKS\n");
}


// void test(){
//   printf("%lf\n", sim_time);
// }

void run()
{
  initialized = 1;
  updateSourceTerm = 1;
  int faceCount = 0;

    for(int ip=0; ip < np; ip++)
    {
      if(demPart[ip].injtime < demTime/timeFactor){
      //if(!demPart[ip].escaped && !demPart[ip].inserted){
        int iIndex = floor((demPart[ip].pos[0]-xmin)/domainDx);
        int jIndex = floor((demPart[ip].pos[1]-ymin)/domainDy);
        int kIndex = floor((demPart[ip].pos[2]-zmin)/domainDz);
        int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;  
        insertToBdBox(ip, cellIndex);//add particle to the cell 
        demPart[ip].inserted = true;
        demPart[ip].active = true;
        rebldNList = true;
      //}
      }
    }
  printf("allocated to BD box\n");
  
  if(rebldNList){
    updateNeighbourList();
    rebldNList = false;
  }

  printf("neighbour list updated\n");
  //rotAng = CURRENT_TIMESTEP*cfdcycles*angVel;
  //printf("timeStep, maxTime %lf %lf\n",timeStep,maxTime);
  demcycles = 0;
  while(demTime < maxTime)
  {
    totalTorque = 0.0; // reset total torque before particle iteration
    demTime += timeStep;
    rotorAngPosition = rotorAngPosition + angVel*timeStep;

      for(int ip=0; ip < np; ip++)
      { 
        //if(demPart[ip].active && !demPart[ip].escaped){
          if(demPart[ip].minWallDist > contactGap*demPart[ip].dia)
          {
            demPart[ip].impactExist = false;
          }
          demPart[ip].minWallDist = 1.0e3; //reset
          forceCalculation(ip);
          maxCnt = fmax(maxCnt,demPart[ip].noOfCnt);
        //}
      }
      //printf("force calculation\n");

      for(int ip=0; ip < np; ip++)
      { 
        //if(demPart[ip].active && !demPart[ip].escaped){
          updatePosition(ip);
        //}
      }

      //printf("stage 2\n");
      for(int ip=0; ip < np; ip++)
      {
          rebldNList = true;
          //reset force for next timestep
          demPart[ip].rotMomc = 0.0;
          demPart[ip].force[0] = 0.0;//demPart[ip].mass;; //gravitational force
          demPart[ip].force[1] = 0.0;//gravitational force
          demPart[ip].force[2] = -demPart[ip].mass; 
          demPart[ip].momentum[0] = 0.0;
          demPart[ip].momentum[1] = 0.0;
          demPart[ip].momentum[2] = 0.0;
      }
   

    if(rebldNList){
      updateNeighbourList();
      rebldNList = false;
    }
    updateSourceTerm = 0;
    //printf("stage 3\n");

    if(ppTotalImpacts > MAXIMPACTS - 100 ||
      totalImpacts > MAXIMPACTS - 100){
      dumpImpact("pw-impact.dat");
      dumpPPImpact("pp-impact.dat");
    }
  
    //}//end if iteration
    demcycles++;
    if(demcycles%ftime == 0)
    {
      demSave();
      writeTorque();
      //dataSave();
      // dumpPPForce();
      // dumpPWForce();
      writeDumpFile("sample");
      //writeCfdCycle("cfdcycle.dat");
      //writeInjectionTime("injtime-");
      dumpImpact("pw-impact.dat");
      dumpPPImpact("pp-impact.dat");
      //printCPUTime();
      printf("Wall time %lf\n",demTime/timeFactor);
    }
  }//end of DEM iteration
    
}

  

