import java.io.*;
class Histogram
{
    static PrintWriter screen = new PrintWriter( System.out, true);// Have to define it, since WriteToDisk uses it for one command
    // these variables have class scope. see Hubbardpage 197 for use of 'protected'
    protected double binsize, binlow, binhigh;
    protected String title;
    protected int SIZE, underflow, overflow,Nsteps;
    
    int[] hist; // define an integer array to store the histogram
    double[] error; //define a double array to store the error on histogram values
    double[] perror; //define a double array to store the percentage error on histogram values

    // contructor method for the class Histogram
    public Histogram(String t, int S, double binlo, double binhi)
    {
        // store the parameters in local variables to be used later
        title = t;
        SIZE = S;
        binlow = binlo;
        binhigh = binhi;        
        //calculate any variables that might be useful later.
        binsize = ( binhigh - binlow)/(double) SIZE;
        hist = new int[SIZE];
        underflow = 0;
        overflow =  0;
        error = new double[SIZE];
        perror = new double[SIZE];

    }
    // ---------------------------------------------------------
    // instance methods start here
    // ---------------------------------------------------------
    public int getSize() { return SIZE;}
    // ---------------------------------------------------------
    public double getBinSize() { return binsize;}
    // ---------------------------------------------------------
    public void fillh( double x)
    {
        if ( x > binlow && x < binhigh)
        {
            // update the correct bin
            int bin = (int)((x - binlow)/binsize);
            hist[bin]++; // add 1 to the bin
            error[bin] =Math.sqrt(hist[bin]);
            perror[bin] = (error[bin]/hist[bin])*100;
        }
        else
        { 
            if (x <= binlow) underflow++;

            if( x>= binhigh) overflow++;

        }

    }
    // -----------------------------------------------------------
    public double getError(int nbin)
    {
        return error[nbin];
    }

    public double getPerror(int nbin)
    {
        return perror[nbin];
    }
    // -----------------------------------------------------------
    public int getUnderflow() {return underflow;}

    public int getOverflow()  {return overflow;}
    // -----------------------------------------------------------
    public String getTitle()
    {
        // return the title of the histogram to the user
        return title;
    }
    // -----------------------------------------------------------
    public int getContent( int nbin)
    { 
        // returns the contnets on bin 'nbin' to the user
        return hist[nbin];
    }

    //-------------------Write to Disk-------------------------
    //Example of a class method to save the histogram data in a file on disk to import 
    //into Excel to make presentation quality graphs and charts
    public void WriteToDisk( double trials) throws IOException
    //public static void WriteToDisk(String t, int s, int[] hist, double[] error, double[] perror, double binlow, double binhigh, double binsize, int underflow, int overflow, double trials, double sumN) throws IOException
    {   //
        //This method handles the writing to disk     
        String filename="..\\"+title+".csv"; //Creates a file with given name in directory one above the class directory
        FileWriter file1 = new FileWriter(filename); //this crates the file 
        PrintWriter outputFile = new PrintWriter (file1); // this sends the output to file1
        // we chose to write the file as a comma separated file (.csv) so you can read it into EXCEL
        screen.println("Writing to disk, please wait....");
        outputFile.println("Title of Histogram: , "+title);
        outputFile.println("Number of trials: , " +trials);
        outputFile.println("Binlow , " +binlow); // note the comma in the text here
        outputFile.println("Binhigh , " +binhigh);
        outputFile.println("Binint , " +binsize);  //ditto the previous comment
        outputFile.println("nbins , " +SIZE); //ditto the previous comment
        outputFile.println("Underflows , "+underflow);
        outputFile.println("Overflows , "+overflow);
        outputFile.println("Bin number , Centre , N , Error , %Error");
        outputFile.println(" , Underflow , "+underflow+" , "+(Math.sqrt(underflow))+" , "+((100*Math.sqrt(underflow))/underflow)); 
        // now make a loop to write the content to each bin to disk, one number at a time
        //together with the x -coordinate of the centre of each bin
        for (int n = 0; n <= SIZE-1; n++)
        {
            //calculate the x coordinates of the centre of each bin
            double binCentre = binlow + binsize/2 + n*binsize;
            outputFile.println(n+" , " + binCentre+" , " +hist[n]+","+error[n]+","+perror[n]);
            //note in the above line we specifically write the comma into the file
        }
        outputFile.println(" , Overflow , "+overflow+" , "+(Math.sqrt(overflow))+" , "+((100*Math.sqrt(overflow))/overflow)); 

        //outputFile.println(" , Sum: ,"+sumN); //Display the sum of all counts in bins
        //outputFile.println(" , Total: ,"+(sumN+underflow+overflow));//Checking that no points were lost, should be equal to a number of trials
        outputFile.close(); // close the output file. THIS IS AN IMPORTANT LINE
        screen.println(" Data written to disk in file " +filename);

    }
    //----------------------Write to disk finished---------------------------
}
