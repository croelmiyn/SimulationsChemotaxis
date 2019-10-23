import ij.IJ;
import ij.gui.*;
import ij.io.*;
import ij.plugin.PlugIn;
import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.nio.ByteBuffer;
import java.util.ArrayList;

public class traj_IO implements PlugIn{

    String dir;
    String fileName;
    String saveName;

    int nFrames;
    int nTraj;
    int nData;

    Traject[] trajectories;
    ArrayList<int[]>[] ptInT;

    boolean canceled;
    String lineSep;

    MersenneTwister rd;
    double D;

    public void run(String arg){

        rd = new MersenneTwister();

        setLoadName();
        readData();

        setSaveName();

        saveFile();

    }


    void setLoadName(){

        OpenDialog od = new OpenDialog("Track_File","");
        dir = od.getDirectory();
        fileName = od.getFileName();

        D = 0.0;

        GenericDialog gd = new GenericDialog("param");
        gd.addNumericField("Measurement_error",D, 1);
        gd.showDialog();

        D = gd.getNextNumber();

    }
    void readData(){

        try {

            FileInputStream file = new FileInputStream(dir+fileName);
            DataInputStream dis = new DataInputStream(file);

            nData = dis.readInt();
            nFrames = dis.readInt();
            nTraj = dis.readInt();

            IJ.log("nb of Traj: "+nTraj);
            IJ.log("nb of Frames: "+nFrames);

            trajectories = new Traject[nTraj];
            for(int t=0; t<nTraj; t++){
                trajectories[t] = new Traject();
            }
            ptInT = new ArrayList[nFrames];
            for(int t=0; t<nFrames; t++){
                ptInT[t] = new ArrayList<int[]>();
            }

            int tn,fn;
            double x,y,phi,pon,meth,cc;
            int[] pointer = new int[2];

            ByteBuffer bb = ByteBuffer.allocate(2*4+6*8);

            for(int i=0; i<nData; i++){
                IJ.showProgress((double)i/nData);

                dis.read(bb.array());

                tn = bb.getInt();
                fn = bb.getInt();
                x = bb.getDouble();
                y = bb.getDouble();
                phi = bb.getDouble();
                pon = bb.getDouble();
                meth = bb.getDouble();
                cc = bb.getDouble();
                bb.clear();

                trajectories[tn].add(fn,x,y,phi,pon,meth,cc);

                pointer = new int[2];
                pointer[0]=tn;
                pointer[1]=trajectories[tn].size()-1;

                ptInT[fn].add(pointer);

            }

        } catch (Exception e){
            IJ.log("Erreur doLoad 1 --> "+e.getMessage());
            IJ.log("Erreur doLoad 2 --> "+e.getCause());
            IJ.log("Erreur doLoad 3 --> "+e.getLocalizedMessage());
        }
        IJ.showStatus("Done");
    }

    private void setSaveName(){

        String sFileName = fileName.replaceFirst(".data",".txt");

        SaveDialog sd = new SaveDialog("Output_File",sFileName,".data");
        saveName = sd.getDirectory()+sd.getFileName();

    }

    private void saveFile(){

        double x,y;
        lineSep = System.lineSeparator();

        String buffer = "traj \t t \t x \t y"+lineSep;

        try {
            FileWriter file = new FileWriter(saveName);
            file.write(buffer);

            for(int i=0; i<nTraj; i++){

                Traject traj = trajectories[i];

                for(int t=0; t<traj.size(); t++){

                    Pt pt = traj.getPt(t);
                    x = pt.getX() + D  * rd.nextGaussian();
                    y = pt.getY() + D  * rd.nextGaussian();

                    buffer = ""+(i+1)+"\t"+pt.getT()+"\t"+x+"\t"+y+lineSep;

                    file.write(buffer);

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
}
