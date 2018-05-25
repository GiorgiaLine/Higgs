import java.io.*;
public class MagnetTrack
{
    static PrintWriter screen = new PrintWriter (System.out, true);
    protected final static double B = 4;
    protected final static double R = 1.2;
    protected static double px, py, angleI, pT, E, pz, Q, Delta, Phi, x, y, z, l, ZAngle;
    
    public static void MagnetTracker( double Energy, double Px, double Py, double Pz, double Charge)
    {
        E = Energy;
        px = Px;
        py = Py;
        pz = Pz;
        Q = Charge;
                
        angleI = Math.atan2(py,px);
        pT = Math.sqrt((px*px)+(py*py));
        ZAngle = Math.atan2(pz,pT);
        Delta = 1000*0.3*B*R*Q/(2*pT);
        Phi = angleI + Delta;
        x = 100*R*Math.cos(Phi);
        y = 100*R*Math.sin(Phi);
        l = Math.sqrt((x*x)+(y*y));
        z = l*Math.tan(ZAngle);
        
        //screen.println("px= " +px+ "py= " +py);
        //screen.println("angleI = " +angleI+ " , pT = " +pT+ ", Delta = " +Delta+ ", Phi = " +Phi+ ", (x,y,z) = (" +x+ "," +y+ "," +z+ ")");
    }
    
    public static double getInitialAngle()
    {
        return angleI;
    }
    
    public static double getInitialMomentum()
    {
        return pT;
    }
    
    public static double getDelta()
    {        
        return Delta;
    }
    
    public static double getPhi()
    {
        //screen.println("phi" + Phi);
        return Phi;
    }
    
    public static double getE()
    {
        return E;
    }
    
    public static double getX()
    {
        return x;
    }
    
    public static double getY()
    {
        return y;
    }
    
    public static double getZAngle()
    {
        return ZAngle;
    }
}
