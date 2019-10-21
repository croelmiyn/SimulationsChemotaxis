import cern.jet.math.Bessel;
import ij.IJ;

public class BesselCalculator {

    private int L;
    private double stp;

    private double[] bk2;
    private double[] bk0;
    private double[] c1;
    private double[] c2;

    BesselCalculator(double maxD, double step){

        stp=step;
        L = (int) (maxD/step) + 1;

        bk2 = new double[L];
        bk0 = new double[L];
        c1 = new double[L];
        c2 = new double[L];
    }

    public void computeCoefs(){

        double z;

        for(int i=0; i<L; i++){

            z=stp*(i+1);

            bk2[i] = Bessel.kn(2,z);
            bk0[i] = Bessel.k0(z);

            c1[i] = z*z/2*(bk0[i]+bk2[i]) - 1;
            c2[i] = 2-z*z*bk2[i];
        }

    }

    public double getC1(double x){

        int i = (int) (x/stp)-1;

        double eps = x/stp-(i+1);

        if(i<0){
            i=0;
            eps=0;
        }
        else if(i>L-2){
            IJ.log(""+x);
            i=L-2;
            eps=0;
        }

        return c1[i] + ((eps!=0)?eps*(c1[i+1]-c1[i]):0);

    }

    public double getC2(double x){

        int i = (int) (x/stp)-1;

        double eps = x/stp-(i+1);

        if(i<0){
            i=0;
            eps=0;
        }
        else if(i>L-2){
            i=L-2;
            eps=0;
        }

        return c2[i] + ((eps!=0)?eps*(c2[i+1]-c2[i]):0);

    }

}
