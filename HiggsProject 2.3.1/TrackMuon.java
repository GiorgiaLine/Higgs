import java.io.*;
import java.util.Random; 
class TrackMuon
/** an accurate 2D tracking of a high energy muon through iron with no magnetic field.*/
{
    // define the identifiers with class scope.
    static BufferedReader keyboard = new BufferedReader (new InputStreamReader(System.in)) ;
    static PrintWriter screen = new PrintWriter( System.out, true);

    static double inputEnergy, stepSize, ironThickness;
    static int NumberOfMuons, Nsteps, numberOfMuons;
    static double px, py, angleI, pT, E, pz, Q, Delta, Phi, x, y, z, l, ZAngle;

    //NumberOfMuons =import3.getNumber();
    
    // Instantiate the class Random, give 'value' class scope
    static Random value = new Random(); 

     static double [] [] Finald1; //static?
     static double [] [] Finald2;
     static double [] [] Finald3;
    
    
    
    static Histogram exitE = new Histogram("Muon exit from iron energy", 50,0,1000);
    static Histogram Detector1 = new Histogram("Detector 1", 50, -205, 25);
    static Histogram Detector2 = new Histogram("Detector 2", 50, -25, 25);
    static Histogram Detector3 = new Histogram("Detector 3", 50, -25, 25);


    //notice defining Histograms here gives the reference pointer (exitE) class scope.
    //--------Class methods start here-------------------------------
    /**   Set the seed and get data from the keyboard */
    private static void getStartingConditions() throws IOException
    {
        // first set the seed for the random number generator so it always
        // produces the same sequence of random numbers
        long seed = 7894694;

        value.setSeed(seed); 
       // screen.println(" Type in starting energy in MeV  ");        
       // inputEnergy = new Double(keyboard.readLine() ).doubleValue();

        screen.println(" Type in a step size in cm  ");
        stepSize = new Double(keyboard.readLine() ).doubleValue();

        ironThickness = 100;

        //screen.println(" Type in the number of muons to track  ");
        numberOfMuons = Import3.getNumber();
        
        NumberOfMuons=numberOfMuons;
        
        
        return;
    }
    //-------------------------------------------------------------------

    /**   Throw a gaussian distribution given the mean and sigma*/
    private static double gauss( double xmean, double sigma )
    {
        // Return a random number with a gaussian distribution
        double newGauss, sum;

        sum=0;
        for (int n=0 ; n<=11; n++)
        {
            sum=sum + value.nextDouble();// use the class Random to make a number
        }
        newGauss = xmean + sigma*(sum -6);

        return newGauss;      

    }
    //------------------------------------------------------------
    /**   Allows output for each muon after the tracking has finished.*/
    private static void lookAtThisMuon(int nsteps, double [][][] track, double finalE, int n)
    {
        // Method to output some information about each muon
        double xlast,ylast,zlast;
        Nsteps=nsteps;
        xlast = track  [nsteps+2][0][n-1];
        ylast = track  [nsteps+2][1][n-1];
        zlast = track  [nsteps+2][2][n-1];

        //screen.println(" energy of muon as it leaves the material = " + finalE +"  MeV");// Have to define finalE from MATHCAD example

        // make histograms  only if the muon left the material
        if (xlast >= ironThickness )
        {
            exitE.fillh(finalE);  // histogram the exit energy
        }
        return;
    }

    //------------ End of class methods --------------------------

    
    public static void TrackerInfo (double DataPrint [] [] []) throws IOException
    {
        //System.out.println(" Event number Num  5, Muon 0, energy" + DataPrint[5][0][0]);
        //System.out.println(" Event number Num  6, Muon 1, Px" + DataPrint[6][1][1]);
        
    }
        
        
    public static void Tracking ( double Data [] [] ) throws IOException
    { 
        // Ask user for input data.
        getStartingConditions();

        //MagnetTrack Test1 = new MagnetTrack( 20000, 26, 1, -6, -1);
        /*
        EnergyLoss iron = new EnergyLoss("iron", 26,55.85,7.87, inputEnergy); //we have to write these classes EnergyLoss & MCS for these constructors to work
        MCS ironMCS = new MCS("iron",26,55.85,7.87, ironThickness, inputEnergy); // Constructor as follows: Name of element, Z, A, density.
        */
        int nsteps; // count number of steps 
        int nmax =200; // maximum allowed number of steps before we stop following a muon.
        // Define a 2D array to store the (x,y) pairs generated as track is followed.
        // allow enough space to store the hit positions on the counters.
        double [] [] [] trackOfMuon = new double[nmax+3] [3] [NumberOfMuons];// see page 169; 
        //nmax +3 to create a slot in array for three detectors when number of steps = nmax         

        double actualMuonEnergy;
        double x; // x coordinate of muon , starts at x=0 just inside iron
        double y; // y coordinate of muon, introduced by multiple scattering
        double z;
        double xnew,ynew,znew;
        double theta; // current angular direction (radians) of muon
        double thetaZ;
        double thetaT;// width of MCS distribution for a muon.
        double Phi;
        double xlast,ylast,zlast,x1,x2,y1,y2,x3,y3,z2,z3,z1, xstart, ystart; 
        
        Finald1 = new double [NumberOfMuons+1] [3];
        Finald2 = new double [NumberOfMuons+1] [3];
        Finald3 = new double [NumberOfMuons+1] [3];
        //double [] [] Finald1 = new double [NumberOfMuons+1] [3];
        //double [] [] Finald2 = new double [NumberOfMuons+1] [3];
        //double [] [] Finald3 = new double [NumberOfMuons+1] [3];

        // Define position and resolution of counters that detect the muon as it 
        // exits the iron.
        double xc1= ironThickness + 10; // x - coord of first counter after the iron
        double xc2 =ironThickness + 20;
        double xc3 =ironThickness + 30;

        double counterYcoordResolution = 250e-5; // sigma of y coord resolution in cms. 
       
        // Start tracking each muon
        //for (int n = 1; n <= NumberOfMuons; n++) 
        for (int n = 1; n <= NumberOfMuons-1; n++) 
        {
            E = Data [n][0]*1e3;
            px = Data [n][1]*1e3;
            py = Data [n][2]*1e3;
            pz = Data [n][3]*1e3;
            Q = Data [n][4];
            //screen.println("\n\n Start tracking muon  : " + E + " , " + px + " , " + py + " , " + pz + " , " + Q );
            
            actualMuonEnergy = E;
            MagnetTrack.MagnetTracker(E,px,py,pz,Q);
            
            //actualMuonEnergy = Data[0][0][0];
            EnergyLoss iron = new EnergyLoss("iron", 26,55.85,7.87, actualMuonEnergy); //we have to write these classes EnergyLoss & MCS for these constructors to work
            MCS ironMCS = new MCS("iron",26,55.85,7.87, ironThickness, actualMuonEnergy, Q); // Constructor as follows: Name of element, Z, A, density.

            x=0; //Set the initial starting position.
            y=0;
            z=0;
            nsteps=0;
            theta = MagnetTrack.getDelta(); // muon starts out parallel to x-axis
            if (theta > Math.PI)
           {
            theta= theta-(Math.PI);
           }
           if (theta < -1*(Math.PI))
           {
            theta= theta+(Math.PI);
           }
            thetaZ = MagnetTrack.getZAngle();
            if (thetaZ > Math.PI)
           {
            thetaZ= thetaZ-(Math.PI);
           }
           if (thetaZ < -1*(Math.PI))
           {
            thetaZ= thetaZ+(Math.PI);
           }
            Phi = MagnetTrack.getPhi();
            if (Phi > Math.PI)
           {
            Phi= Phi-(Math.PI);
           }
           if (Phi < -1*(Math.PI))
           {
            Phi= Phi+(Math.PI);
           }
            //screen.println("\n\n Start tracking muon  " + n + " ,energy =  " + actualMuonEnergy );
            // In this program we are working in units of cms.

            while ( x < ironThickness && nsteps < nmax ) // Note the 2 conditions here
            {
                // Step is the direction in the x-direction. If the muon is scattered by and angle
                // theta then the amount of material the muon travels through is d = step/cos(theta)                
                double step = Math.min( stepSize, ironThickness-x);
                // Ensure the final step just reaches the end of the iron, crucial line.               

                // Find width of MCS distribution for this muon travelling a distance stepSize 
                // through material.
                //thetaT= ironMCS.mcsTheta0(actualMuonEnergy, step);
                thetaT= ironMCS.getThetaL(actualMuonEnergy, step);

                // Generate a random  angle with mean 0 with gaussian spread to add to current 
                // direction.
                //screen.println("theta before=" +theta+ ", "+thetaT+ ", "+actualMuonEnergy+ ", "+step);
                theta= theta + gauss (0,thetaT);
                //screen.println("theta=" +theta);
                thetaZ= thetaZ + gauss (0,thetaT);
                
                double dxy = step/Math.cos(theta);
                double d = dxy/Math.cos(thetaZ);
                //double d = step/Math.cos(theta);

                // Find energy loss going through d cm of material.
                actualMuonEnergy = actualMuonEnergy - iron.getEnergyL(actualMuonEnergy) *  d; 
                // Warning: the above line assumes that the energy loss can be regarded as being
                // essentially constant for the muon travelling a distance  'step'.
                // If this is not true then it is necessary to change step.

                if(actualMuonEnergy < 0)
                {
                    screen.print(" Energy of muon goes negative.. abandon it");
                    break; // This causes the 'for' loop to terminate.
                }

                xnew = x + step ; // calculate next (x,y) position.
                ynew = y + d*Math.sin(theta); 
                znew = z + d*Math.sin(thetaZ);
                //screen.println(" Number of Muons "+NumberOfMuons);

                //screen.println(" tracking.. nsteps " + nsteps + " xnew " + xnew + " ynew: " + ynew + " znew: " + znew);
                screen.flush();
                //String anykey;
                //anykey = keyboard.readLine();// pause until any key is pressed.
                // Store these co-ordinates           
                
                trackOfMuon  [nsteps] [0][n-1] = xnew;
                trackOfMuon  [nsteps] [1][n-1] = ynew;
                trackOfMuon  [nsteps] [2][n-1] = znew;

                // Update coordinates in order to take the next step.
                x = xnew;
                y = ynew;
                z = znew;

                nsteps++;
                // At this point will return to the start of the 'while' loop and take another step.

                if(nsteps == nmax) screen.println(" Too many steps for muon " + n + ",  abandon it");

            }
            // Finished tracking this muon, do some analysis on the results, and calculate hit
            // coordinates on the counters
            //screen.println("angles : "+theta+", "+Phi);
            double tanAngle = Math.tan(theta)+Math.tan(Phi);
            double AngleFinal = Math.atan(tanAngle);
             if (AngleFinal > Math.PI)
           {
            AngleFinal= AngleFinal-(Math.PI);
           }
           if (AngleFinal < -1*(Math.PI))
           {
            AngleFinal= AngleFinal+(Math.PI);
           }
            double ThetaFinalZ = 0;//thetaZ;
            //screen.println("angles combined : "+AngleFinal);
            
            //screen.println("angle1 " +tanAngle+ " , angle2 " + AngleFinal+ " , angle3 " +ThetaFinalZ);
            
            xstart = MagnetTrack.getX();
            ystart = MagnetTrack.getY();
            
            double xFinal = x + xstart; // final coords on track
            double yFinal = y + ystart;
            double zFinal = z;

            // Work out y-hit position on each counter and SMEAR it by the resolution.
            double yhitOnC1 = (xc1 - xFinal)*Math.tan(AngleFinal) + yFinal;
            yhitOnC1=gauss(yhitOnC1,counterYcoordResolution);
            double zhitOnC1 = (xc1 - xFinal)*Math.tan(ThetaFinalZ) + zFinal;
            zhitOnC1=gauss(zhitOnC1,counterYcoordResolution);

            double yhitOnC2= (xc2 - xFinal)*Math.tan(AngleFinal) +yFinal; 
            yhitOnC2 = gauss( yhitOnC2,counterYcoordResolution);
            double zhitOnC2= (xc2 - xFinal)*Math.tan(ThetaFinalZ) + zFinal;
            zhitOnC2 = gauss( zhitOnC2,counterYcoordResolution);

            double yhitOnC3= (xc3 - xFinal)*Math.tan(AngleFinal) +yFinal; 
            yhitOnC3 = gauss( yhitOnC3,counterYcoordResolution);
            double zhitOnC3= (xc3 - xFinal)*Math.tan(ThetaFinalZ) + zFinal;
            zhitOnC3 = gauss( zhitOnC3,counterYcoordResolution);
            
            double xreal1= x*Math.cos(AngleFinal) + 10; // x - coord of first counter after the iron
            double xreal2 =x*Math.cos(AngleFinal) + 20;
            double xreal3 =x*Math.cos(AngleFinal) + 30;
            
            // Add these coords into the array;      n-1 because we want to account for the first entry in the array
            // n-1 , because array starts from 0
            trackOfMuon  [nsteps] [0][n-1] = xreal1;
            trackOfMuon  [nsteps] [1][n-1] = yhitOnC1;
            trackOfMuon  [nsteps] [2][n-1] = zhitOnC1;
            //Detector1.fillh(yhitOnC1);
            trackOfMuon  [nsteps +1] [0][n-1] = xreal2;
            trackOfMuon  [nsteps +1] [1][n-1] = yhitOnC2;
            trackOfMuon  [nsteps +1] [2][n-1] = zhitOnC2;
            //Detector2.fillh(yhitOnC2);
            trackOfMuon  [nsteps +2] [0][n-1] = xreal3;
            trackOfMuon  [nsteps +2] [1][n-1] = yhitOnC3;
            trackOfMuon  [nsteps +2] [2][n-1] = zhitOnC3;
            //Detector3.fillh(yhitOnC3);
            // pass the data to this method for any further processing
            lookAtThisMuon(nsteps,trackOfMuon,actualMuonEnergy, n); //
               
            xlast = trackOfMuon  [nsteps+2][0][n-1];
            ylast = trackOfMuon  [nsteps+2][1][n-1];
            zlast = trackOfMuon  [nsteps+2][2][n-1];
        
            x2 = trackOfMuon  [nsteps+1][0][n-1];
            y2 = trackOfMuon  [nsteps+1][1][n-1];
            z2 = trackOfMuon  [nsteps+1][2][n-1];
        
            x1 = trackOfMuon  [nsteps][0][n-1];
            y1 = trackOfMuon  [nsteps][1][n-1];
            z1 = trackOfMuon  [nsteps][2][n-1];

            Finald3 [n-1][0] = xlast;
            Finald3 [n-1][1] = ylast;
            Finald3 [n-1][2] = zlast;
        
            Finald2 [n-1][0] = x2;
            Finald2 [n-1][1] = y2;
            Finald2 [n-1][2] = z2;
        
            Finald1 [n-1][0] = x1;
            Finald1 [n-1][1] = y1;
            Finald1 [n-1][2] = z1;
        
            //screen.println(" last (x,y,z) of track =  ( " + xlast + " , " + ylast + " , " +zlast +")" );
            // Now generate the next muon
            //}
        }
        //WritePositionToDisk(Nsteps+2, trackOfMuon);  //for 1 muon leave inside for loop, outside loo for more than 1 muon
        // -----------------------------------------------------------
         /*void Linefit(double [][] xcoord,double [][] ycoord,int n,float* Ans) {
             // -----------------------------------------------------------
             if(n < 3) {
                Ans[0] = 0.0;
                Ans[1] = 0.0;
                return;
            }
            float Count = 0.0;
            float Sumx  = 0.0;
            float Sumy  = 0.0;
            float Sumxy = 0.0;
            float Sumxx = 0.0;
            float Sumyy = 0.0;
            for(int j=1; j<=n; j++) {
                if(y[j-1] != 0.0) {
                  Sumx = Sumx + x[j-1];
                  Sumy = Sumy + y[j-1];
                  Count= Count+ 1.0;
                }
            }
            if(Count <= 1.0) {
                       Ans[0] = 0.0;
                       Ans[1] = 0.0;
                       return;
                    }
                    float Ymed = Sumy/Count;
                    float Xmed = Sumx/Count;
                    for(int j=1; j<=n; j++) {
                        if(y[j-1] != 0.0) {
                            float Scartx = x[j-1] - Xmed;
                            float Scarty = y[j-1] - Ymed;
                            Sumxy  = Sumxy + Scartx*Scarty;
                            Sumxx  = Sumxx + Scartx*Scartx;
                   Sumyy  = Sumyy + Scarty*Scarty;
                }
            }
            // -----------------------------------------------
            // Fit Parameters:
            // -----------------------------------------------
            if(Sumxx == 0.0) {
                       Ans[0] = 0.0;
                       Ans[1] = 0.0;
                       return;
                    }
                    float A = Sumxy/Sumxx;
                    float B = Ymed - A*Xmed;
                    float E = 0.0;
                    if(Count >= 3.0) {
                       E = (Sumyy - Sumxy*A)/(Count-2.0);
                    }
                    Ans[0] = A;
                    Ans[1] = B;
                }
                // --------------------------------------------------------------------
                float xxx[3];
                float yyy[3];
                float Ans[2];
                // --------------------------------------------------------------------
                // Test Linefit:
                // --------------------------------------------------------------------
                // xxx[0] = 1.0;
                // yyy[0] = 4.0*xxx[0] + 3.5;
                // xxx[1] = 2.0;
                // yyy[1] = 4.0*xxx[1] + 3.5;
                // xxx[2] = 3.0;
                // yyy[2] = 4.0*xxx[2] + 3.5;
                // int n3 = 3;
                // Linefit(xxx,yyy,n3,Ans);
                // cout<<" Gradient "<<Ans[0]<<endl;
                // cout<<" Intercept "<<Ans[1]<<endl;*/

        //screen.println( "x1 " + Finald1[0][0] + "  y1 " + Finald1[0][1] + "  z1 " +Finald1[0][2]);
        //screen.println( "x2 " + Finald2[0][0] + "  y2 " + Finald2[0][1] + "  z3 " +Finald2[0][2]);
        //screen.println( "x3 " + Finald3[3][0] + "  y3 " + Finald3[4][1] + "  z3 " +Finald3[3][2]);
        //Detector1.WriteToDisk(NumberOfMuons);
        //Detector2.WriteToDisk(NumberOfMuons);
        //Detector3.WriteToDisk(NumberOfMuons);
        
        exitE.WriteToDisk(NumberOfMuons);           

        //GetFinald1();   
        //GetFinald2(); 
        //GetFinald3(); 
        // All muons done, finish program. If necessary, write histograms to disk at this point
    }
    //--------------Writting track of mouns into .csv file---------------------------------------------
   

    public static double [] [] GetFinald1()
    {
        return Finald1; 
    }
    public static double [] [] GetFinald2()
    {
        return Finald2; 
    }
    public static double [] [] GetFinald3()
    {
        return Finald3; 
    }
    
    
    
    public static void WritePositionToDisk(double nsteps, double [][][]trackOfMuon) throws IOException
    {
        String filename="..\\Track.csv"; //Creates a file with given name in directory one above the class directory
        FileWriter file1 = new FileWriter(filename); //this crates the file 
        PrintWriter outputFile = new PrintWriter (file1);
        
        outputFile.println("Moun number , X position , Y position, Z position");
        //Double for loop to access all data in 3D array
        for (int j = 1; j<=NumberOfMuons; j++)
        {
            for (int i = 0; i<=nsteps; i++)
            {
                //screen.println(j + " || " + trackOfMuon[i][0][j-1] + " || " + trackOfMuon[i][1][j-1]); // For debugging
                outputFile.println( j + "," + trackOfMuon[i][0][j-1] + "," + trackOfMuon[i][1][j-1] + "," +trackOfMuon[i][2][j-1]);
            }
        }        
        outputFile.close();
    }
    
    /*public static void Arrays(double nmuons, double [][]Detector) throws IOException
    {
        //String filename="..\\Track.csv"; //Creates a file with given name in directory one above the class directory
        //FileWriter file1 = new FileWriter(filename); //this crates the file 
        //PrintWriter outputFile = new PrintWriter (file1);
        
        outputFile.println("Moun number , X position , Y position, Z position");
        //Double for loop to access all data in 3D array
        for (int j = 1; j<=NumberOfMuons; j++)
        {
            for (int i = 0; i<=nsteps; i++)
            {
                //screen.println(j + " || " + trackOfMuon[i][0][j-1] + " || " + trackOfMuon[i][1][j-1]); // For debugging
                outputFile.println( j + "," + trackOfMuon[i][0][j-1] + "," + trackOfMuon[i][1][j-1] + "," +trackOfMuon[i][2][j-1]);
            }
        }        
        outputFile.close();
    }*/
    //--------------Finished Writting track-------------------------------------------------------------
}