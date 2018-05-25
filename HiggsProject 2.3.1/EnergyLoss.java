import java.io.*;
class EnergyLoss
{
    static PrintWriter screen=new PrintWriter (System.out, true);
    
    protected String material;
    protected double Z, A, Ro, E, b, p, y, I, energyLoss, EL;
    protected final double mu= 106;
    protected final double me= 0.511;
    
    public double getEnergyL( double E)
    { 
          
        p = Math.sqrt(E*E-mu*mu);
        b = (p/(Math.sqrt((p*p)+(mu*mu))));
        y = (1/(Math.sqrt(1-(b*b))));
        I= (1.35e-5)*Z;
        
        energyLoss = 0.307*Ro*(Z/A)*(1/(b*b))*((Math.log((2*me*b*b*y*y)/I))-(b*b));
        return energyLoss;
    }
    
     public EnergyLoss(String material, double z, double a, double ro, double e)
    {
        
        Z=z;
        A=a;
        Ro=ro;
        E=e;
        /*        
        double b = (p/(Math.sqrt((p*p)+(mu*mu))));
        double y = (1/(Math.sqrt(1-(b*b))));
        double I= (1.35e-5)*z;
        */
        //double energyLoss = 0.307*ro*(z/a)*(1/(b*b))*((Math.log((2*me*b*b*y*y)/I))-(b*b));
        EL=getEnergyL(e);
        //screen.println(" Energy loss "+EL);
        //screen.println(" momentum "+p);
        //screen.println(" beta "+b);
        //screen.println(" gamma "+y);
        
    }
        
    public double getBeta()
    {
        return b;
    }
    
}