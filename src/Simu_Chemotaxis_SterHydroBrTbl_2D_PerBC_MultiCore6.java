//////////////////////////////////////////////////////////////////////////////////////////////////
// ImageJ Plugin for simulation of chemotaxis with steric and hydrodynamic interactions   		//
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
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import mpi.rc.IJ.IJutilities.ConcurrencyUtils;
import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class Simu_Chemotaxis_SterHydroBrTbl_2D_PerBC_MultiCore6 implements PlugIn{

    /* Main Plugin for simulation of chemotaxis of interacting rods with hydrodynamic interactions (2D motion + force dipole + Hele-Shaw approximation)
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
    private double h;                   // channel height
    private double dip;                 // hydrodynamic dipole factor
    private double lf;


    private double kappa;               // normalized elasticity
    private double kin;               	// normalized inelastic coeff
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

    private double zrot;                  // accessory param = 3/L^2

    // for computations
    boolean runs[];

    double PhiM[];
    double dxM[];
    double dyM[];
    double vP[];

    double vxt[];
    double vyt[];
    double wt[];

    double uxM[];
    double uyM[];


    MersenneTwister rd2 = new MersenneTwister();

    /* neighboring relations
    ArrayList<ArrayList<Integer>> neighbors = new ArrayList<ArrayList<Integer>>(); // lists of neighboring pixels, within rad, for each pixel
    ArrayList<ArrayList<Integer>> nbgParts = new ArrayList<ArrayList<Integer>>(); // lists of cells in the neighborhood, within rad, for each pixel
    //*/
    // concentration:

    double grad;
    double avgCC;


    ImagePlus imp;

    ImageProcessor ip;

    StackWindow iw;

    String lineSep;

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
        kappa = 10;
        kin = 0.1;
        sTS = 100;

        h = 20;
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
        gd.addNumericField("E/eta*l^3/2 : ", kappa, 0);
        gd.addNumericField("Inelastic K /eta : ", kin, 2);
        gd.addNumericField("Length_of_Flag (px): ", lf, 0);
        gd.addNumericField("Rot_diff_coeff (fr-1): ", dr, 0);
        gd.addNumericField("Channel_height (microns): ", h, 0);
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
        kappa = (double) gd.getNextNumber();
        kin = (double) gd.getNextNumber();
        lf = (double) gd.getNextNumber();
        dr = (double) gd.getNextNumber();
        h = (double) gd.getNextNumber();
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

        kappa /= sTS*Math.PI*L/Math.log(L/e); // dt E ln(L/e)/(pi eta L l^3/2)
        kin /=sTS*Math.PI*L/Math.log(L/e); // dt Kin ln(L/e)/(pi eta L )

        zrot = 3/(L*L); // for rotational friction

        Ltot = L + lf;   // the flagellum length is variable
        dip = h*L/(6*Math.log(L/e)); // Ltot/2 is dipole length
        sqrtK = Math.sqrt(12)/h;
        //BesselCalculator bc = new BesselCalculator(frameSizeX,0.01);
        //dip = L*L/4 / ( bc.getC1(sqrtK*L/2) + bc.getC2(sqrtK*L/2) ); // the flow speed is the swimming speed at the from tip of the cell

        //dr = 8*Math.log(Ltot/(e*L/Ltot+0.1*8/Ltot))/(Ltot*Ltot*Ltot)/100/sTS; //rotational diffusion (2Dr): 100 corresponds to 1 frame = 0.01s
        dr *= 2/sTS;
        dr = Math.sqrt(dr); // sqrt(2Dr)

        dTheta = Math.sqrt(dTheta/sTS);

        Ltot/=2; //for computations: half length

        Network.setTimeStep(1.0);
        Network.setDissConstants(Ka_off,Ka_on,Ks_off,Ks_on);
        Network.setAdaptPrecision(1.0);
        Network.setRelativeCheA(1.0);
        Network.setRelativeCheRCheB(1.0,1.0);
        Network.setRelativeCheYCheZ(1.0,1.0);
        Network.setNreceptors(nTar,nTsr);

        Network.setIniMethState(computeSteadyMethylation(avgCC)); // 0 to 8

        setSaveName();

        nPairs = (nObjects*(nObjects-1))/2;
        pairs = new int[nPairs][2];
        int p=0;
        for(int i=0; i<nObjects; i++){
            for(int j=i+1; j<nObjects; j++){
                pairs[p][0] = i;
                pairs[p][1] = j;
                p++;
            }
        }

    }

    private void Compute(){

        // parallel computing utils
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        int p = nObjects / nthreads;
        if(p<1){
            p=1;
            nthreads = nObjects;
        }
        int p2 = nPairs / nthreads;
        Future<?>[] futures = new Future[nthreads];

        locks = new Lock[nObjects];
        for(int k=0;k<nObjects;k++){
            locks[k] = new ReentrantLock();
        }

        BC = new BesselCalculator[nthreads];
        for (int k=0;k<nthreads;k++) {
            BC[k] = new BesselCalculator(frameSizeX,0.01);
            BC[k].computeCoefs();
        }


        // main data
        posX = new double[nObjects][nFrames];
        posY = new double[nObjects][nFrames];
        phi = new double[nObjects][nFrames];

        vP = new double[nObjects];   // static rods velocity = 0
        runs = new boolean[nObjects];

        long dt1=0;
        long dt2=0;
        long dt3=0;
        long t0,t1;

        // Initialize the neighboring relations

        //hydro: not relevant

        // initialize the positions

        MersenneTwister rd1 = new MersenneTwister();

        // for increment computations
        PhiM = new double[nObjects];
        dxM = new double[nObjects];
        dyM = new double[nObjects];
        uxM = new double[nObjects];
        uyM = new double[nObjects];

        vxt = new double[nObjects];
        vyt = new double[nObjects];
        wt = new double[nObjects];

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


        // make the rods parallel

        // parallelization block
        for (int l = 0; l < nthreads; l++) {

            final int firstK = l * p2;
            final int lastK = ( l == (nthreads-1))? (nPairs) : (firstK + p2);

            futures[l] = ConcurrencyUtils.submit(new Runnable()
            {
                public void run(){

                    for(int k=firstK; k<lastK; k++){

                        // insert here the code

                        computeMechForces(pairs[k],0);
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
        } catch (ExecutionException ex) {
            IJ.log("Execution Exception");
            IJ.log("Message --> "+ex.getMessage());
            IJ.log("Cause --> "+ex.getCause());
            IJ.log("LocMessage --> "+ex.getLocalizedMessage());
        }
        //end parallelisation block

        // Position implementation needs to be separate
        // parallelization block
        for (int l = 0; l < nthreads; l++) {

            final int firstK = l * p;
            final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);

            futures[l] = ConcurrencyUtils.submit(new Runnable()
            {
                public void run(){

                    for(int k=firstK; k<lastK; k++){

                        // insert here the code

                        updatePositionAndOrientationInit(k);

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
        } catch (ExecutionException ex) {
            IJ.log("Execution Exception");
            IJ.log("Message --> "+ex.getMessage());
            IJ.log("Cause --> "+ex.getCause());
            IJ.log("LocMessage --> "+ex.getLocalizedMessage());
        }
        //end parallelisation block


        // initialize velocity, becomes dynamic
        //rd1 = new Random();

        for(int k=0; k<nObjects; k++){

            vP[k] = v0 + dv0 * rd1.nextGaussian();//Math.sqrt(-2*Math.log(rd1.nextDouble())) * Math.cos(2 * Math.PI * rd1.nextDouble());
            vP[k]/=sTS;
        }

        // initialization of the chemotactic pathway
        runs = new boolean[nObjects];
        chemNetwork = new Network[nObjects];

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

            //rd2 = new Random();

            IJ.showProgress((double)t/nFrames);

            for(int k=0; k<nObjects; k++){

                phi[k][t] = phi[k][t-1];

                posX[k][t] = posX[k][t-1] ;
                posY[k][t] = posY[k][t-1] ;


            }

            final int tL = t; // for par computing

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
                IJ.log("Chem and run");
            } catch (ExecutionException ex) {
                IJ.log("Execution Exception");
                IJ.log("Message --> "+ex.getMessage());
                IJ.log("Cause --> "+ex.getCause());
                IJ.log("LocMessage --> "+ex.getLocalizedMessage());
                IJ.log("Chem and run");
            }
            //end parallelisation block

            //* parallelization block
            for (int l = 0; l < nthreads; l++) {

                final int firstK = l * p2;
                final int lastK = ( l == (nthreads-1))? (nPairs) : (firstK + p2);
                final int thread =l;

                futures[l] = ConcurrencyUtils.submit(new Runnable()
                {
                    public void run(){

                        for(int k=firstK; k<lastK; k++){

                            // insert here the code

                            computeHydroForces(pairs[k],tL, thread);
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
                IJ.log("Forces");
            } catch (ExecutionException ex) {
                IJ.log("Execution Exception");
                IJ.log("Message --> "+ex.getMessage());
                IJ.log("Cause --> "+ex.getCause());
                IJ.log("LocMessage --> "+ex.getLocalizedMessage());
                IJ.log("Forces");
            }
            //end parallelisation block */

            t1 = System.currentTimeMillis();
            dt1 += t1-t0;

            for(int kk=0; kk<sTS; kk++){

                for (int l = 0; l < nthreads; l++) {

                    final int firstK = l * p;
                    final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);

                    futures[l] = ConcurrencyUtils.submit(new Runnable()
                    {
                        public void run(){

                            for(int k=firstK; k<lastK; k++){

                                incrementRunHydro(k,tL);

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

                //* parallelization block
                for (int l = 0; l < nthreads; l++) {

                    final int firstK = l * p2;
                    final int lastK = ( l == (nthreads-1))? (nPairs) : (firstK + p2);

                    futures[l] = ConcurrencyUtils.submit(new Runnable()
                    {
                        public void run(){

                            for(int k=firstK; k<lastK; k++){

                                // insert here the code

                                computeMechForces(pairs[k],tL);
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
                    IJ.log("Forces");
                } catch (ExecutionException ex) {
                    IJ.log("Execution Exception");
                    IJ.log("Message --> "+ex.getMessage());
                    IJ.log("Cause --> "+ex.getCause());
                    IJ.log("LocMessage --> "+ex.getLocalizedMessage());
                    IJ.log("Forces");
                }
                //end parallelisation block */

                // Position implementation needs to be separate
                // parallelization block
                for (int l = 0; l < nthreads; l++) {

                    final int firstK = l * p;
                    final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);

                    futures[l] = ConcurrencyUtils.submit(new Runnable()
                    {
                        public void run(){

                            for(int k=firstK; k<lastK; k++){

                                // insert here the code

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
            }

            // reset the cells in the neighborhoods

            //*/
            t0 = System.currentTimeMillis();
            dt2 += t0-t1;

            // saving
            Traject traj;
            for(int k=0; k<nObjects; k++){
                traj = trajectories.remove(trajNb[k]);
                traj.add(t,posX[k][t],posY[k][t],phi[k][t],chemNetwork[k].P_on,chemNetwork[k].meth,conc(k,t));
                trajectories.add(trajNb[k],traj);
            }


        }

        IJ.log("trajectory calculated...");
        IJ.log("part chem+Hydro ="+((double)dt1/(dt1+dt2+dt3)) );
        IJ.log("part updatePos+steric ="+((double)dt2/(dt1+dt2+dt3)) );
        IJ.log("total time ="+((double)(dt1+dt2+dt3)/1000) );
    }

    private void generateMovie(){

        // generates the film

        imp = IJ.createImage("Simu_Chemotaxis_HydroSteric_"+parm,"32-bit", frameSizeX, frameSizeY, nFrames);

        ip = imp.getProcessor();

        double px, dr2, cz,X,Y, ck,sk;
        int x0,y0;

        double a = 1.0;

        //int[][] black = new int[frameSizeX][frameSizeY];
        float[][] intens = new float[frameSizeX][frameSizeY];

        for(int t=0; t<nFrames; t++){

            imp.setSlice(t+1);
            //ip.setIntArray(black);

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

                        intens[modulo((x0+xx),frameSizeX)][modulo((yy+y0),frameSizeY)]+= 255 * Math.exp(-a*dr2) * cz;;

                        //ip.putPixel((x0+xx)%frameSizeX,(yy+y0)%frameSizeY,(int)px);
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

                    /*
                    dos.writeInt(i);
                    dos.writeInt(pt.getT());
                    dos.writeDouble(pt.getX());
                    dos.writeDouble(pt.getY());
                    dos.writeDouble(pt.getPhi());
                    dos.writeDouble(pt.getPon());
                    dos.writeDouble(pt.getMeth());
                    dos.writeDouble(pt.getC());
                    //*/
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

        String buffer = "Parameters Hydro Steric"+lineSep;

        buffer +="frameSizeX\t"+frameSizeX+lineSep;
        buffer +="frameSizeY\t"+frameSizeY+lineSep;
        buffer +="nFrames\t"+nFrames+lineSep;
        buffer +="nObjects\t"+nObjects+lineSep;
        buffer +="v0\t"+v0+lineSep;
        buffer +="dv0\t"+dv0+lineSep;
        buffer +="L\t"+L+lineSep;
        buffer +="e\t"+e+lineSep;
        buffer +="kappa\t"+kappa+lineSep;
        buffer +="kin\t"+kin+lineSep;
        buffer +="lf\t"+lf+lineSep;
        buffer +="dr\t"+dr+lineSep;
        buffer +="h\t"+h+lineSep;
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

    private double sgn(double x){
        if(x>0.0){
            return 1.0;
        }
        else if(x<0.0){
            return -1.0;
        }
        else{
            return 0.0;
        }
    }

    private double abs(double x){
        return Math.abs(x);
    }

    private double force(double x){
        return kappa * Math.pow(x,1.5); // 4x avt  Math.min(kappa * Math.pow(x,1.5),e)
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
        double bias = Math.pow(chemNetwork[k].cheYp,hillCoef); // run[t] = f(chemPathway[t-1])
        chemNetwork[k].updateMWCmodel(conc(k,t));

        if(runs[k]){

            // do we stay in the run phase in the next step ?
            double prt = Math.exp(-rtR*bias);
            synchronized(rd2){runs[k]=(rd2.nextDouble()<prt);}
        }
        else{

            // do we stay tumble in the next step ?
            double ptr = Math.exp(-trR);// /bias . changed 06 01 2016
            synchronized(rd2){runs[k]=!(rd2.nextDouble()<ptr);}  // since rd2.nextDouble()<ptr : proba to stay tumble.

        }

        //prepare displacement computation

        uxM[k] = Math.cos(phi[k][t]);
        uyM[k] = Math.sin(phi[k][t]);

        wt[k] =0;

        vxt[k] = (runs[k])?(vP[k] * uxM[k]):0;
        vyt[k] = (runs[k])?(vP[k] * uyM[k]):0;

    }

    void computeHydroForces(int[] index, int t, int thread){

        int k=index[0];
        int cl = index[1];

        double cosPhii = uxM[k];
        double sinPhii = uyM[k];

        double dfxi = -cosPhii*Ltot;
        double dfyi = -sinPhii*Ltot;

        double cosPhij,sinPhij, dx,dy, dfxj,dfyj;

        double r2;

        double dXi,dYi,dXj,dYj,dPhii,dPhij;
        double[] Us;

        dx=posX[cl][t] - posX[k][t];
        if(abs(dx)>frameSizeX/2){ // periodic boundary conditions
            dx -= sgn(dx)*frameSizeX;
        }
        dy=posY[cl][t] - posY[k][t];
        if(abs(dy)>frameSizeY/2){ // periodic boundary conditions
            dy -= sgn(dy)*frameSizeY;
        }

        cosPhij = uxM[cl];
        dfxj = -cosPhij*Ltot;
        sinPhij = uyM[cl];
        dfyj = -sinPhij*Ltot;

        r2 = dx*dx+dy*dy;

        // hydrodynamics
        if(r2>e*e) { //separate

            Us = getSpeeds(-dx,-dy,cosPhii,sinPhii,cosPhij,sinPhij,vP[k],vP[cl],runs[k],runs[cl],thread); //uxi,uyi,uxj,uyj
            dXi = Us[0];
            dYi = Us[1];
            dXj = Us[2];
            dYj = Us[3];

            dPhii = cosPhii * Us[1] - sinPhii * Us[0];
            dPhij = cosPhij * Us[3] - sinPhij * Us[2];

            Us = getSpeeds(-dx+dfxi,-dy+dfyi,cosPhii,sinPhii,cosPhij,sinPhij,-vP[k],vP[cl],runs[k],runs[cl],thread); //uxi,uyi,uxj,uyj
            dXi += Us[0];
            dYi += Us[1];
            dXj += Us[2];
            dYj += Us[3];

            dPhii -= cosPhii * Us[1] - sinPhii * Us[0];
            dPhij += cosPhij * Us[3] - sinPhij * Us[2];

            Us = getSpeeds(-dx-dfxj,-dy-dfyj,cosPhii,sinPhii,cosPhij,sinPhij,vP[k],-vP[cl],runs[k],runs[cl],thread); //uxi,uyi,uxj,uyj
            dXi += Us[0];
            dYi += Us[1];
            dXj += Us[2];
            dYj += Us[3];

            dPhii += cosPhii * Us[1] - sinPhii * Us[0];
            dPhij -= cosPhij * Us[3] - sinPhij * Us[2];

            Us = getSpeeds(-dx+dfxi-dfxj,-dy+dfyi-dfyj,cosPhii,sinPhii,cosPhij,sinPhij,-vP[k],-vP[cl],runs[k],runs[cl],thread); //uxi,uyi,uxj,uyj
            dXi += Us[0];
            dYi += Us[1];
            dXj += Us[2];
            dYj += Us[3];

            dPhii -= cosPhii * Us[1] - sinPhii * Us[0];
            dPhij -= cosPhij * Us[3] - sinPhij * Us[2];

            dPhii/=Ltot;
            dPhij/=Ltot;

            dXi /= 2;
            dYi /= 2;
            dXj /= 2;
            dYj /= 2;

            safeAddV(k, dXi, dYi);
            safeAddV(cl, dXj, dYj);
            safeAddW(k, dPhii);
            safeAddW(cl, dPhij);


        }

    }

    double[] getSpeeds(double dx,double dy,double ci,double si,double cj,double sj,double vi, double vj, boolean ri, boolean rj, int thread){

        double r=Math.sqrt(dx*dx+dy*dy);
        double cor =dip / (r*r);

        double c1 = BC[thread].getC1(sqrtK*r);
        double c2 = BC[thread].getC2(sqrtK*r);

        double[] out = new double[4];

        double ux=dx/r;
        double uy=dy/r;

        double ct, st, vx, vy;
        if(rj) {
            /*
            ct = cj * ux + sj * uy;
            st = cj * uy - sj * ux;

            vx = c1 + c2 * ct * ct;
            vy = c2 * ct * st;

            out[0] = vj * (vx * cj - vy * sj) * cor;
            out[1] = vj * (vx * sj + vy * cj) * cor;
            //*/
            ct = cj * ux + sj * uy;
            ct*=c2;

            out[0] = vj * cor * (c1 * cj + ct * ux);
            out[1] = vj * cor * (c1 * sj + ct * uy);

        }
        if(ri) {
            /*
            ct = ci * (-ux) + si * (-uy);
            st = ci * (-uy) - si * (-ux);

            vx = c1 + c2 * ct * ct;
            vy = c2 * ct * st;

            out[2] = vi * (vx * ci - vy * si) * cor;
            out[3] = vi * (vx * si + vy * ci) * cor;
            */

            ct = ci * ux + si * uy;
            ct*=c2;

            out[2] = vi * cor * (c1 * ci + ct * ux);
            out[3] = vi * cor * (c1 * si + ct * uy);
        }
        return out;

    }

    private void incrementRunHydro(int k, int t){

        dxM[k] = vxt[k] ;
        dyM[k] = vyt[k] ;

        PhiM[k] = wt[k];
        synchronized(rd2){PhiM[k] += dr * rd2.nextGaussian();} // rotational diffusion

        if(!runs[k]){

            // reorientation in the tumble phase
            synchronized(rd2){PhiM[k] += dTheta *(rd2.nextDouble()-0.5);} //* Math.sqrt(-2*Math.log(rd.nextDouble())) * Math.cos(2 * Math.PI * rd.nextDouble())

        }

    }

    void computeMechForces(int[] index,int t){

        int k=index[0];
        int cl = index[1];

        double dx,dy;

        dx = posX[cl][t] - posX[k][t];
        if (abs(dx) > frameSizeX / 2) { // periodic boundary conditions
            dx -= sgn(dx) * frameSizeX;
        }
        dy = posY[cl][t] - posY[k][t];
        if (abs(dy) > frameSizeY / 2) { // periodic boundary conditions
            dy -= sgn(dy) * frameSizeY;
        }

        // repulsion part
        if (dx * dx + dy * dy < L * L) {

            double cosPhii = uxM[k];
            double sinPhii = uyM[k];
            double cosPhij = uxM[cl];
            double sinPhij = uyM[cl];

            double[] out = getStericRepulsion(dx, dy, cosPhii, sinPhii, cosPhij, sinPhij, vP[k], vP[cl]);

            safeAddPos(k, out[0], out[1]);
            safeAddPos(cl, out[2], out[3]);
            safeAddOrr(k, out[4]);
            safeAddOrr(cl, out[5]);

        } // if r2


    }

    private double[] getStericRepulsion(double dx,double dy,double cosPhii, double sinPhii,double cosPhij,double sinPhij, double vi,double vj){

        double cji, drui,druj, aa, bb, vx,vy, fij, epsi1, epsj1, d, vui, dX,dY,dPhi;

        double[] out = new double[6]; //dxk,dyk,dxcl,dycl,dphik,dphicl  k=i, cl=j

        cji = cosPhii * cosPhij + sinPhii * sinPhij; // cos(phi_j-phi_i)

        // Contact Check

        if (cji * cji == 1.0) { // parallel rods

            drui = dx * cosPhii + dy * sinPhii;
            druj = -dx * sinPhii + dx * cosPhii; // misused = drvi in fact

            if (abs(drui) < (L) && abs(druj) < (e)) { // contact

                d = e - abs(druj);
                fij = force(d); // saturation of the response
                fij /= 4; // Zperp*fij

                // elastic repulsion force
                aa = fij * sgn(druj) * sinPhii;
                bb = fij * sgn(druj) * cosPhii;
                // inelastic friction force
                aa += kin * (cosPhij * vj - cosPhii * vi);
                bb -= kin * (sinPhij * vj - sinPhii * vi);

                //safeAddPos(k, aa, -bb);
                out[0] += aa;
                out[1] -= bb;

                // neighbor implementation
                //safeAddPos(cl, -aa, bb);
                out[2] -= aa;
                out[3] += bb;


            }

        } else { // non parallel rods

            drui = dx * cosPhii + dy * sinPhii;
            druj = dx * cosPhij + dy * sinPhij;

            // intersection des droites
            epsi1 = (drui - cji * druj) / (1 - cji * cji);
            epsj1 = (-druj + cji * drui) / (1 - cji * cji);

            if (abs(epsi1) <= L / 2 && abs(epsj1) <= L / 2) {// crossing rods

                // weird alignment correction

                d = Math.sqrt(2.0-2.0*cji);
                if(cji>0.0) {
                    vx = (cosPhii + cosPhij) / d;
                    vy = (sinPhii + sinPhij) / d;
                }else{
                    vy = (cosPhii + cosPhij) / d;
                    vx = -(sinPhii + sinPhij) / d;
                }

                //*        Force computation
                if ((dx * vx + dy * vy) > 0) { // dr.v must be < 0
                    vx = -vx;
                    vy = -vy;
                }
                //fij = force(e); // saturation of the response
                fij = e / 2;

                // Inelastic friction force Fin = kin ( vP(j) uj.vPerp - vP(i) ui.vPerp ) vPerp
                aa = kin * (vj * (-cosPhij * vy + sinPhij * vx) - vi * (-cosPhii * vy + sinPhii * vx));
                bb = aa * (vx); // component y
                aa *= (-vy);  // component x

                // total force
                vx *= fij;
                vx += aa;
                vy *= fij;
                vy += bb; // v<-(fij*v+Fin)

                //* dPhi = Zrot * (epsi ui) x (+ fij v + Fin)
                dPhi = zrot * epsi1 * (cosPhii * vy - sinPhii * vx);
                out[4]+=dPhi;
                //*/
                //* dPhj = Zrot * (epsj uj) x (- fij v - Fin)
                dPhi = zrot * epsj1 * (-cosPhij * vy + sinPhij * vx);
                out[5]+=dPhi;
                //*/


                // part impl
                vx *= 0.25;
                vy *= 0.25; // v<- 0.25*v
                //dri+=  ( (Zpar-ZPerp) *   (fij*v+Fin).ui  *  ui +  Zperp * (fij*v+Fin)), with (Zpar-ZPerp)=Zperp = 0.25
                vui = (vx * cosPhii + vy * sinPhii);
                dX = (vui * cosPhii + vx);
                dY = (vui * sinPhii + vy);
                //safeAddPos(k, dX, dY);
                out[0]+=dX;
                out[1]+=dY;

                // neighbor impl:
                //drj+= fij * ( (Zpar-ZPerp) *   (-(fij*v+Fin)).uj      *      uj   +  Zperp * (-(fij*v+Fin))), with (Zpar-ZPerp)=Zperp
                vui = (vx * cosPhij + vy * sinPhij);
                dX = -(vui * cosPhij + vx);
                dY = -(vui * sinPhij + vy);
                //safeAddPos(cl, dX, dY);
                out[2]+=dX;
                out[3]+=dY;



            /* Simple realignment
            PhiM[k] += d;
            //*/
            } else { // non parallel, non crossing rods
                // compute minimal distance points actually on the rod
                if (abs(epsi1) >= abs(epsj1)) {
                    epsi1 = sgn(epsi1) * L / 2;
                    epsj1 = epsi1 * cji - druj;
                    if (abs(epsj1) > L / 2) {
                        epsj1 = sgn(epsj1) * L / 2;
                    }
                }
                if (abs(epsj1) >= abs(epsi1)) {
                    epsj1 = sgn(epsj1) * L / 2;
                    epsi1 = epsj1 * cji + drui; // + since dr is antisymmetric by i-j exchange
                    if (abs(epsi1) > L / 2) {
                        epsi1 = sgn(epsi1) * L / 2;
                    }
                }

                // compute vector between min distance points
                vx = dx + epsj1 * cosPhij - epsi1 * cosPhii;  // v in the direction opposit to the force
                vy = dy + epsj1 * sinPhij - epsi1 * sinPhii;

                d = Math.sqrt(vx * vx + vy * vy);

                if (d < e) { // contact without crossing

                    vx /= d;
                    vy /= d;

                    fij = force(e - d); // saturation of the response

                    // Inelastic friction force Fin = kin ( vP(j) uj.vPerp - vP(i) ui.vPerp ) vPerp
                    aa = kin * (vj * (-cosPhij * vy + sinPhij * vx) - vi * (-cosPhii * vy + sinPhii * vx));
                    bb = aa * (vx); // component y
                    aa *= (-vy);  // component x

                    // total force
                    vx *= fij;
                    vx -= aa;
                    vy *= fij;
                    vy -= bb; // v<-(fij*v-Fin)

                    // dPhi = Zrot * (epsi ui) x (-fij v + Fin)
                    dPhi = -zrot * epsi1 * (cosPhii * vy - sinPhii * vx);
                    //safeAddOrr(k, dPhi);
                    out[4]+=dPhi;

                    // dPhi = Zrot * (epsj uj) x (+fij v - Fin)
                    dPhi = zrot * epsj1 * (cosPhij * vy - sinPhij * vx);
                    //safeAddOrr(cl, dPhi);
                    out[5]+=dPhi;

                    vx *= 0.25;
                    vy *= 0.25;

                    //dri +=  ( (Zpar-ZPerp) *   (-fij *v+Fin).ui      *      ui   +  Zperp * (-fij*v+Fin)), with (Zpar-ZPerp)=Zperp = 0.25
                    vui = (vx * cosPhii + vy * sinPhii);
                    dX = -(vui * cosPhii + vx);
                    dY = -(vui * sinPhii + vy);
                    //safeAddPos(k, dX, dY);
                    out[0]+=dX;
                    out[1]+=dY;

                    // neighbor impl
                    //drj += fij * ( (Zpar-ZPerp) *   (+fij *v-Fin).uj      *      uj   +  Zperp * (+fij*v-Fin))
                    vui = (vx * cosPhij + vy * sinPhij);
                    dX = (vui * cosPhij + vx);
                    dY = (vui * sinPhij + vy);
                    //safeAddPos(cl, dX, dY);
                    out[2]+=dX;
                    out[3]+=dY;


                }

            }

        }// else non parallel

        return out;

    }

    private void updatePositionAndOrientation(int k, int t){

        // position and orientation iteration
        phi[k][t] = phi[k][t] + PhiM[k];

        posX[k][t] = posX[k][t] + dxM[k] ;
        posY[k][t] = posY[k][t] + dyM[k] ;

        // update of the orientation and speed
        if(runs[k]) {
            vxt[k] -=  (vP[k] * uxM[k]);
            vyt[k] -=  (vP[k] * uyM[k]);
        }

        uxM[k] = Math.cos(phi[k][t]);
        uyM[k] = Math.sin(phi[k][t]);

        if(runs[k]) {
            vxt[k] +=  (vP[k] * uxM[k]);
            vyt[k] +=  (vP[k] * uyM[k]);
        }

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

    private void updatePositionAndOrientationInit(int k){
        int t=0;
        // position and orientation iteration
        phi[k][t] = phi[k][t] + PhiM[k];

        posX[k][t] = posX[k][t] + dxM[k] ;
        posY[k][t] = posY[k][t] + dyM[k] ;

        // crossing boundary conditions
        if(posX[k][t]>frameSizeX){
            posX[k][t] = posX[k][t]-frameSizeX;
        }
        if(posX[k][t]<0){
            posX[k][t] = frameSizeX+posX[k][t];
        }

        //  periodic boundary conditions Y
        if(posY[k][t]>frameSizeY){
            posY[k][t] = posY[k][t]-frameSizeY;
        }
        if(posY[k][t]<0){
            posY[k][t] = frameSizeY+posY[k][t];
        }

    }

    private void safeAddPos(int k, double dX, double dY){
        locks[k].lock();
        try{
            dxM[k] += dX;
            dyM[k] += dY;
        }
        finally{
            locks[k].unlock();
        }
    }
    private void safeAddOrr(int k, double dPhi){
        locks[k].lock();
        try{
            PhiM[k] += dPhi;
        }
        finally{
            locks[k].unlock();
        }
    }

    private void safeAddV(int k, double dX, double dY){
        locks[k].lock();
        try{
            vxt[k] += dX;
            vyt[k] += dY;
        }
        finally{
            locks[k].unlock();
        }
    }
    private void safeAddW(int k, double dPhi){
        locks[k].lock();
        try{
            wt[k] += dPhi;
        }
        finally{
            locks[k].unlock();
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