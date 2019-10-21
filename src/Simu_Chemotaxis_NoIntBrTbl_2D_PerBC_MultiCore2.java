//////////////////////////////////////////////////////////////////////////////////////////////////
// ImageJ Plugin for simulation of chemotaxis without interactions								//
// Uses Concurrency Utilities extracted of Parallel colt by Piotr Wendykier 					//
// Uses the Mersenne Twister algorithm by (c) Michael Lecuyer for java (v13)					//
// CC BY SA	by Remy Colin and the Max Planck Society      										//
//////////////////////////////////////////////////////////////////////////////////////////////////
// Date   : 2019-10-21															//
// Author : Rémy Colin															//
// Email  : remycolin@laposte.net / remy.colin@synmikro.mpi-marburg.mpg.de		//
// Adress : Max Planck Institute for Terrestrial Microbiology					//
// 			Karl-von-Frisch-strasse 10											//
//          35043 Marburg, Germany												//
//////////////////////////////////////////////////////////////////////////////////
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.StackWindow;
import ij.io.SaveDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import mpi.rc.IJ.IJutilities.ConcurrencyUtils;
import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class Simu_Chemotaxis_NoIntBrTbl_2D_PerBC_MultiCore2 implements PlugIn{

    /* Rods without interactions
     */

    private double posX[][], posY[][], phi[][];  // pour mémoriser la position des différents objets sur les différentes slices
    private int nObjects;               // le nombre d'objets sur lesquels faire le calcul
    private int nFrames;
    private int frameSizeX;
    private int frameSizeY;
    private double dr; 					//sqrt(2 Dr)
    private double dv0;
    private double v0;
    private double L;					// length of the particles
    private double e;                   // half-thickness of the particles
    private double Ltot;                // cell + flagellum
    private double lf;

    private double sTS;					// separation of Compression vs velocity time scales

    private double rtR;
    private double trR;
    private double dTheta;

    private double Ka_on;
    private double Ka_off;
    private double Ks_on;
    private double Ks_off;
    private double hillCoef;
    private int nTar;
    private int nTsr;
    private double adaptRate;

    Network[] chemNetwork;

    // for computations
    boolean runs[];

    double PhiM[];
    double dxM[];
    double dyM[];
    double vP[];

    double uxM[];
    double uyM[];


    MersenneTwister rd2 = new MersenneTwister();

    private int coarseSizeX,coarseSizeY;

    double grad;
    double avgCC;


    ImagePlus imp;

    ImageProcessor ip;

    StackWindow iw;

    String lineSep;

    ResultsTable rt;

    private String parm;

    ArrayList<Traject> trajectories;
    int[] trajNb;

    String saveName,saveNameLog;

    Lock[] locks; // for safe access to arrays

    BesselCalculator[] BC;
    double sqrtK;

    int[][] pairs;
    int nPairs;


    public void run(String arg0) {

        Initialise();

        Compute();

        // save the trajectories

        saveLog();

        saveData();

        IJ.log("true trajectory saved...");

        generateMovie();

        IJ.log("done");

    }// run

    private void Initialise(){

        lineSep = System.getProperty("line.separator");

        nFrames = 5000;
        nObjects = 1000;               // le nombre d'objets sur lesquels faire le calcul
        frameSizeX = 256;
        frameSizeY = 256;

        v0 = 0.2;
        dv0=0.02;

        L=4;
        e=1;

        sTS = 100;

        dr = 0.001;
        lf = 8;


        rtR = 0.01;
        trR = 0.10;
        dTheta = 0.12;
        Ka_on = 3000;
        Ka_off = 20; // MeAsp
        Ks_on = 1e4;
        Ks_off = 1e3;
        hillCoef = 10;
        nTar = 10;
        nTsr = 10;
        adaptRate = 0.01;

        grad = 0.00025;
        avgCC= 100;

        parm = ""+System.currentTimeMillis();

        GenericDialog gd = new GenericDialog("nom du fichier de sortie");
        gd.addNumericField("BoxSize_X_direction : ", frameSizeX, 0);
        gd.addNumericField("BoxSize_Y_direction : ", frameSizeY, 0);
        gd.addNumericField("Length_of_the_Movie_(frames) : ", nFrames, 0);
        gd.addNumericField("Number_of_Objects : ", nObjects, 0);
        gd.addNumericField("Norm_of_the_velocity : ", v0, 2);
        gd.addNumericField("Std_Dev_of_the_norm_of_the_velocity : ", dv0, 2);
        gd.addMessage("Elastic Interaction Parameters");
        gd.addNumericField("Length_of_Rods (px) : ", L, 0);
        gd.addNumericField("Thickness_of_Rods (px) : ", e, 0);
        gd.addNumericField("Length_of_Flag (px): ", lf, 0);
        gd.addNumericField("Rot_diff_coeff (fr-1): ", dr, 0);
        gd.addNumericField("Time_Scales_Sep : ", sTS, 0);
        gd.addMessage("Chemotaxis Pathway Parameters");
        gd.addNumericField("Run_to_Tumble_rate_(fr^-1) : ", rtR, 2);
        gd.addNumericField("Tumble_to_Run_rate_(fr^-1) : ", trR, 2);
        gd.addNumericField("Tumbling Rotational rate (rad/fr) :", dTheta, 2);
        gd.addNumericField("Koff_tar (uM):", Ka_off, 2);
        gd.addNumericField("Kon_tar (uM):", Ka_on, 2);
        gd.addNumericField("Koff_tsr (uM):", Ks_off, 2);
        gd.addNumericField("Kon_tsr (uM):", Ks_on, 2);
        gd.addNumericField("Hill_coef_motor :", hillCoef, 1);
        gd.addNumericField("nTar :", nTar, 0);
        gd.addNumericField("nTsr :", nTsr, 0);
        gd.addNumericField("adapt_rate (fr-1):", adaptRate, 2);
        gd.addMessage("Chemical gradient Parameters");
        gd.addNumericField("gradSlope (px-1):", grad, 5);
        gd.addNumericField("avgC (uM):", avgCC, 0);
        gd.addMessage("Save option");
        gd.addStringField("label:", parm, 0);
        gd.showDialog();
        frameSizeX = (int) gd.getNextNumber();
        frameSizeY = (int) gd.getNextNumber();
        nFrames = (int) gd.getNextNumber();
        nObjects = (int) gd.getNextNumber();
        v0 = (double) gd.getNextNumber();
        dv0 = (double) gd.getNextNumber();
        L = (double) gd.getNextNumber();
        e = (double) gd.getNextNumber();
        lf = (double) gd.getNextNumber();
        dr = (double) gd.getNextNumber();
        sTS = (double) gd.getNextNumber();

        rtR = (double) gd.getNextNumber();
        trR = (double) gd.getNextNumber();
        dTheta = (double) gd.getNextNumber();
        Ka_off = (double) gd.getNextNumber();
        Ka_on = (double) gd.getNextNumber();
        Ks_off = (double) gd.getNextNumber();
        Ks_on = (double) gd.getNextNumber();
        hillCoef = (double) gd.getNextNumber();
        nTar = (int) gd.getNextNumber();
        nTsr = (int) gd.getNextNumber();
        adaptRate = (double) gd.getNextNumber();

        grad = (double) gd.getNextNumber();
        avgCC = (double) gd.getNextNumber();

        parm = gd.getNextString();

        Ltot = L + lf;   

        //dr = 8*Math.log(Ltot/(e*L/Ltot+0.1*8/Ltot))/(Ltot*Ltot*Ltot)/100/sTS; //rotational diffusion (2Dr): 100 corresponds to 1 frame = 0.01s
        dr/= sTS;
        dr = Math.sqrt(2*dr); // sqrt(2Dr)

        Ltot/=2; //for computations: half length

        Network.setTimeStep(1.0/sTS);
        Network.setDissConstants(Ka_off,Ka_on,Ks_off,Ks_on);
        Network.setAdaptPrecision(1.0);
        Network.setRelativeCheA(1.0);
        Network.setRelativeCheRCheB(1.0,1.0);
        Network.setRelativeCheYCheZ(1.0,1.0);
        Network.setNreceptors(nTar,nTsr);

        Network.setIniMethState(computeSteadyMethylation(avgCC)); // 0 to 8

        setSaveName();



    }

    private void Compute(){

        // parallel computing utils
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        int p = nObjects / nthreads;
        if(p<1){
            p=1;
            nthreads = nObjects;
        }
        Future<?>[] futures = new Future[nthreads];

        locks = new Lock[nObjects];
        for(int k=0;k<nObjects;k++){
            locks[k] = new ReentrantLock();
        }




        // main data
        posX = new double[nObjects][nFrames];
        posY = new double[nObjects][nFrames];
        phi = new double[nObjects][nFrames];

        vP = new double[nObjects];   // static rods velocity = 0

        long dt1=0;
        long dt2=0;
        long dt3=0;
        long t0,t1;

        // Initialize the neighboring relations

        // hydro: not relevant

        // initialize the positions

        MersenneTwister rd1 = new MersenneTwister();

        // for increment computations
        PhiM = new double[nObjects];
        dxM = new double[nObjects];
        dyM = new double[nObjects];
        uxM = new double[nObjects];
        uyM = new double[nObjects];

        for(int k=0; k<nObjects; k++){

            // initial position at random
            posX[k][0] = rd1.nextDouble()* frameSizeX;
            posY[k][0] = rd1.nextDouble()* frameSizeY;
            phi[k][0] = rd1.nextDouble()* 2*Math.PI;
            uxM[k] = Math.cos(phi[k][0]);
            uyM[k] = Math.sin(phi[k][0]);

            // periodic boundary conditions
            posX[k][0] -= frameSizeX * Math.floor(posX[k][0]/frameSizeX);
            posY[k][0] -= frameSizeY * Math.floor(posY[k][0]/frameSizeY);

        }



        // initialize velocity, becomes dynamic
        //rd1 = new Random();

        for(int k=0; k<nObjects; k++){

            vP[k] = v0 + dv0 * rd1.nextGaussian();//Math.sqrt(-2*Math.log(rd1.nextDouble())) * Math.cos(2 * Math.PI * rd1.nextDouble());
            vP[k] /= sTS;

        }

        // initialization of the chemotactic pathway
        runs = new boolean[nObjects];
        chemNetwork = new Network[nObjects];

        double dt = 1.0/sTS;

        rtR*=dt;
        trR*=dt;
        dTheta*=dt;
        dTheta = Math.sqrt(dTheta);

        for(int k=0; k<nObjects; k++){

            // initial state at random
            equilibrateNetwork(k,0);

        }

        // trajectory initiation for save

        trajectories = new ArrayList<Traject>();
        trajNb = new int[nObjects];
        for(int k=0; k<nObjects; k++){
            Traject traj = new Traject();
            traj.add(0,posX[k][0],posY[k][0],phi[k][0],chemNetwork[k].P_on,chemNetwork[k].meth,conc(k,0));
            trajectories.add(traj);
            trajNb[k] = trajectories.size()-1;
        }

        IJ.log("initialize done..");

        for(int t=1; t<nFrames; t++){

            IJ.showProgress((double)t/nFrames);

            for(int k=0; k<nObjects; k++){

                phi[k][t] = phi[k][t-1];

                posX[k][t] = posX[k][t-1] ;
                posY[k][t] = posY[k][t-1] ;


            }

            final int tL = t; // for parallel computing

            for(int kk=0; kk<sTS; kk++){ // compression time scale

                t0 = System.currentTimeMillis();

                // parallelization block
                for (int l = 0; l < nthreads; l++) {

                    final int firstK = l * p;
                    final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);

                    futures[l] = ConcurrencyUtils.submit(new Runnable()
                    {
                        public void run(){

                            for(int k=firstK; k<lastK; k++){

                                // insert here the code

                                advanceChemAndRun(k,tL);

                                updatePositionAndOrientation(k,tL);

                                // end of calculus
                            }

                        }
                    });
                }
                try {
                    ConcurrencyUtils.waitForCompletion(futures);
                } catch (InterruptedException ex) {
                    IJ.log("Interruption Exception");
                    IJ.log("Message --> "+ex.getMessage());
                    IJ.log("Cause --> "+ex.getCause());
                    IJ.log("LocMessage --> "+ex.getLocalizedMessage());
                    IJ.log("Update Or Pos");
                } catch (ExecutionException ex) {
                    IJ.log("Execution Exception");
                    IJ.log("Message --> "+ex.getMessage());
                    IJ.log("Cause --> "+ex.getCause());
                    IJ.log("LocMessage --> "+ex.getLocalizedMessage());
                    IJ.log("Update Or Pos");
                }
                //end parallelisation block

                t1 = System.currentTimeMillis();
                dt3 += t1-t0;


            }// for kk Compression

            // saving
            Traject traj;
            for(int k=0; k<nObjects; k++){
                traj = trajectories.remove(trajNb[k]);
                traj.add(t,posX[k][t],posY[k][t],phi[k][t],chemNetwork[k].P_on,chemNetwork[k].meth,conc(k,t));
                trajectories.add(trajNb[k],traj);
            }


        }

        IJ.log("trajectory calculated...");
        IJ.log("part chem+forces ="+((double)dt1/(dt1+dt2+dt3)) );
        IJ.log("part reinit neighbors ="+((double)dt2/(dt1+dt2+dt3)) );
        IJ.log("part updatePosAndNeighbors ="+((double)dt3/(dt1+dt2+dt3)) );
        IJ.log("total time ="+((double)(dt1+dt2+dt3)/1000) );
    }

    private void generateMovie(){

        // generates the film

        imp = IJ.createImage("Simu_Chemotaxis_HydroSteric_"+parm,"32-bit", frameSizeX, frameSizeY, nFrames);

        ip = imp.getProcessor();

        double px, dr2, cz,X,Y, ck,sk;
        int x0,y0;

        double a = 1.0;
		
        float[][] intens = new float[frameSizeX][frameSizeY];

        for(int t=0; t<nFrames; t++){

            imp.setSlice(t+1);
            
            intens = new float[frameSizeX][frameSizeY];

            for(int k=0; k<nObjects; k++){

                x0 = (int) posX[k][t];
                y0 = (int) posY[k][t];

                cz = 0.25;

                ck = Math.cos(phi[k][t]);
                sk = Math.sin(phi[k][t]);

                for(int xx=-(int)(L+e); xx<=(L+e); xx++){

                    for(int yy=-(int)(L+e); yy<=(L+e); yy++){

                        X = (xx+x0-posX[k][t]) * ck + (yy+y0-posY[k][t]) * sk;
                        Y = -(xx+x0-posX[k][t]) * sk + (yy+y0-posY[k][t]) * ck;

                        dr2 = 4*X*X/L/L + 4*Y*Y/e/e;

                        intens[modulo((x0+xx),frameSizeX)][modulo((yy+y0),frameSizeY)]+= 255 * Math.exp(-a*dr2) * cz;
						
                    }

                }

            }
            ip.setFloatArray(intens);

            IJ.showProgress((double) t /nFrames);

        }
        imp.setSlice(1);

        iw = new StackWindow(imp);

        IJ.log("done");


    }


    private void setSaveName(){

        String sFileName = "Trajectories_SimuDensityChemotaxis_"+parm+"";;

        SaveDialog sd = new SaveDialog("Output_File",sFileName,".data");
        saveName = sd.getDirectory()+sd.getFileName();

        sFileName = "Log_SimuDensityChemotaxis_"+parm+"";

        sd = new SaveDialog("Log_File","",sFileName,".txt");
        saveNameLog = sd.getDirectory() + sd.getFileName();


    }
    private void saveData(){

        Traject traj;
        Pt pt;

        int nTraj = trajectories.size();
        int nData = 0;
        for(int i=0; i<nTraj; i++){
            nData+=trajectories.get(i).size();
        }

        try {
            FileOutputStream file = new FileOutputStream(saveName);

            // standard format = nb of data arrays, nFrames, nObjects, t[], data1[][], data2[][], ...

            DataOutputStream dos = new DataOutputStream(file);
            ByteBuffer bb = ByteBuffer.allocate(2*4+6*8);
            //bb.order(ByteOrder.nativeOrder());

            dos.writeInt(nData);
            dos.writeInt(nFrames);
            dos.writeInt(nTraj);

            for(int i=0; i<nTraj; i++){
                IJ.showProgress((double)i/nTraj);
                traj = trajectories.get(i);
                for(int j=0; j<traj.size(); j++){

                    pt = traj.getPt(j);

                    bb.putInt(i);
                    bb.putInt(pt.getT());
                    bb.putDouble(pt.getX());
                    bb.putDouble(pt.getY());
                    bb.putDouble(pt.getPhi());
                    bb.putDouble(pt.getPon());
                    bb.putDouble(pt.getMeth());
                    bb.putDouble(pt.getC());

                    dos.write(bb.array());

                    bb.clear();

                    
                }
            }
            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");

    }

    private void saveLog(){

        String buffer = "Parameters No Int rot Diff"+lineSep;

        buffer +="frameSizeX\t"+frameSizeX+lineSep;
        buffer +="frameSizeY\t"+frameSizeY+lineSep;
        buffer +="nFrames\t"+nFrames+lineSep;
        buffer +="nObjects\t"+nObjects+lineSep;
        buffer +="v0\t"+v0+lineSep;
        buffer +="dv0\t"+dv0+lineSep;
        buffer +="L\t"+L+lineSep;
        buffer +="e\t"+e+lineSep;
        buffer +="lf\t"+lf+lineSep;
        buffer +="dr\t"+dr+lineSep;
        buffer +="sTS\t"+sTS+lineSep;

        buffer +="rtR\t"+rtR+lineSep;
        buffer +="trR\t"+trR+lineSep;
        buffer +="dTheta\t"+dTheta+lineSep;
        buffer +="Ka_off\t"+Ka_off+lineSep;
        buffer +="Ka_on\t"+Ka_on+lineSep;
        buffer +="Ks_off\t"+Ks_off+lineSep;
        buffer +="Ks_on\t"+Ks_on+lineSep;
        buffer +="hillCoef\t"+hillCoef+lineSep;
        buffer +="nTar\t"+nTar+lineSep;
        buffer +="nTsr\t"+nTsr+lineSep;
        buffer +="adaptRate\t"+adaptRate+lineSep;

        buffer +="grad\t"+grad+lineSep;
        buffer +="avgCC\t"+avgCC+lineSep;

        try {
            FileWriter file = new FileWriter(saveNameLog);
            file.write(buffer);

            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");

    }


    private int modulo(int x,int n){
        int a=x%n;
        return (a<0)?(a+n):a;
    }

    private double conc(int k, int t){

        return avgCC * (1 + grad * (2*posX[k][t] - frameSizeX)); // cmax/2 + cmax*slope*(x-xm) = avgc (1+2*slope*(x-xm))

    }

    private void equilibrateNetwork(int k, int t){

        double bias, ptr,prt;
        MersenneTwister rd1 = new MersenneTwister();

        Network.setIniMethState(computeSteadyMethylation(conc(k,t)));
        chemNetwork[k] = new Network();
        chemNetwork[k].setAdaptRate(adaptRate);
        chemNetwork[k].updateMWCmodel(conc(k,t));
        bias = Math.pow(chemNetwork[k].cheYp,hillCoef);

        for(int t1=0; t1<100; t1++){

            if(runs[k]){

                prt = Math.exp(-rtR*bias);

                // do we stay in the runs phase in the next step ?
                runs[k]=(rd1.nextDouble()<prt);
            }
            else{

                ptr = Math.exp(-trR/bias);

                // do we stay tumble in the next step ?

                runs[k]=!(rd1.nextDouble()<ptr);  // since rd2.nextDouble()<ptr : proba to stay tumble.
            }
        }

    }

    private double computeSteadyMethylation(double cc){

        double eps_m = (Math.log(2) - nTar*Math.log( (1+cc/Ka_off)/(1+cc/Ka_on) ) - nTsr * Math.log((1+cc/Ks_off)/(1+cc/Ks_on))) / (nTar+nTsr);
        double initM;

        if(eps_m>1.0) initM=0; // Endres & Wingreen, 2006, piece-wise linear
        else if(eps_m>0.0) initM= (1.0-eps_m)/0.5;
        else if (eps_m>-0.6) initM = 2.0-eps_m/0.3;
        else if  (eps_m>-1.1) initM = 4.0-(eps_m+0.6)/0.25;
        else if (eps_m>-2.0) initM = 6.0-(1.1+eps_m)/0.9;
        else if (eps_m>-3.0) initM = 7.0-(2.0+eps_m);
        else initM=8.0;

        return initM;
    }

    private void advanceChemAndRun(int k, int t){

        // update network
        chemNetwork[k].updateMWCmodel(conc(k,t));
        double bias = Math.pow(chemNetwork[k].cheYp,hillCoef);

        //prepare displacement computation
        dxM[k] = 0;
        dyM[k] = 0;

        PhiM[k]=0;
        synchronized(rd2){PhiM[k]= dr * rd2.nextGaussian();} // rotational diffusion

        uxM[k] = Math.cos(phi[k][t]);
        uyM[k] = Math.sin(phi[k][t]);

        if(runs[k]){

            // straight movement in the run phase
            dxM[k] += vP[k] * uxM[k];
            dyM[k] += vP[k] * uyM[k];

            // do we stay in the run phase in the next step ?
            double prt = Math.exp(-rtR*bias);
            synchronized(rd2){runs[k]=(rd2.nextDouble()<prt);}
        }
        else{

            // reorientation in the tumble phase
            synchronized(rd2){PhiM[k]+= dTheta *(rd2.nextDouble()-0.5);} //* Math.sqrt(-2*Math.log(rd.nextDouble())) * Math.cos(2 * Math.PI * rd.nextDouble())

            // do we stay tumble in the next step ?
            double ptr = Math.exp(-trR);// /bias 
            synchronized(rd2){runs[k]=!(rd2.nextDouble()<ptr);}  // since rd2.nextDouble()<ptr : proba to stay tumble.

        }

    }

    private void updatePositionAndOrientation(int k, int t){

        // position and orientation iteration
        phi[k][t] = phi[k][t] + PhiM[k];

        posX[k][t] = posX[k][t] + dxM[k] ;
        posY[k][t] = posY[k][t] + dyM[k] ;

        // crossing boundary conditions
        if(posX[k][t]>frameSizeX){
            posX[k][t] = posX[k][t]-frameSizeX;
            beginNewTraj(k,t,true);

        }
        if(posX[k][t]<0){
            posX[k][t] = frameSizeX+posX[k][t];
            beginNewTraj(k,t,true);

        }

        //  periodic boundary conditions Y
        if(posY[k][t]>frameSizeY){
            posY[k][t] = posY[k][t]-frameSizeY;
            beginNewTraj(k,t,false);
        }
        if(posY[k][t]<0){
            posY[k][t] = frameSizeY+posY[k][t];
            beginNewTraj(k,t,false);
        }

    }


    private synchronized void beginNewTraj(int k, int t, boolean refreshPathway){

        if(refreshPathway){
            Network.setIniMethState(computeSteadyMethylation(conc(k,t)));
            chemNetwork[k] = new Network();
            chemNetwork[k].setAdaptRate(adaptRate);
        }

        Traject traj =  trajectories.get(trajNb[k]);
        if(!traj.isEmpty()){
            traj = new Traject();
            trajectories.add(traj);
            trajNb[k] = trajectories.size()-1;
        }
    }


}