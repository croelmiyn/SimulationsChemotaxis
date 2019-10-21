import ij.IJ;
import ij.io.OpenDialog;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;

public class DataFileLoader {

    int nData,nFrames,nTraj;

    public String dir,fileName;

    ArrayList<int[]>[] ptInT;
    Traject[] trajectories;
    boolean canceled;

    public DataFileLoader(){
        canceled = false;
        setLoadName();
    }
    public DataFileLoader(String dir,String fileName){
        canceled = false;
        this.dir = dir;
        this.fileName = fileName;
    }

    public void setLoadName(){
        OpenDialog od = new OpenDialog("Track_File","");
        dir = od.getDirectory();
        fileName = od.getFileName();
        canceled = (fileName==null);
    }
    public void readData(){

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


    public Traject[] getTrajectories(){
        return trajectories;
    }

    public ArrayList<int[]>[] getPtInT() {
        return ptInT;
    }

    public boolean wasCanceled(){
        return canceled;
    }
}
