import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.measure.Calibration;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Calculate_Intensities implements PlugIn 
{
    
    double radius=5;
    double radius_angstrom=0.7;
    double scale=-1;//no scale
    int imheight;
    double fov=-1;
    Point[] points;
    ImagePlus imp;
    RoiManager roim;
    ArrayList<Double> meanarray= new ArrayList<Double>();
    
    public void run(String arg) 
    {
            
            imp = IJ.getImage(); 
            imheight=imp.getHeight();
            Roi roi=imp.getRoi();
            if (roi==null || roi.getType()!=Roi.POINT)
            {
                IJ.log("Please create only a point selection");
                return;
            }
            points = roi.getContainedPoints();
            int len=points.length;
            
            Calibration cal = imp.getCalibration();
            String unit = cal.getUnit();
            scale = getScale(unit);
            if (scale!=-1)
            {
                fov=cal.getX(imheight);
                radius=radius_angstrom*imheight/(fov*scale);
            }
            
            if(DoDialog())
            {
                roim=RoiManager.getRoiManager();
                //new RoiManager();
                int[] xarray= new int[len];
                int[] yarray= new int[len];
                for (int ind=0; ind<len;++ind)
                {
                    xarray[ind]=points[ind].x-(int)radius/2;
                    yarray[ind]=points[ind].y-(int)radius/2;
                }
                
                for (int ind=0; ind<len;++ind)
                {
                    //OvalRoi oroi=new OvalRoi(xarray[0],yarray[0]-radius/2,radius,radius);
                    OvalRoi oroi=optimizeRoi(xarray[ind],yarray[ind]);
                    //imp.setRoi(oroi);
                    
                }
            }
            double mean=calculateMean(meanarray);
            double std=calculateStandardDeviation(meanarray);
            IJ.log("Number of atoms found:\t"+len);
            IJ.log("Mean of atom intensities:\n"+mean);
            IJ.log("Standard deviation of means:\n"+std);
            imp.setRoi(roi);

    }

    public boolean DoDialog()
    {
        boolean baccepted=false;
        while (!baccepted)
        {
            NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Radius");
            gd.setSmartRecording(true);
            gd.addMessage("choose radius for roi...");
            if (scale!=-1)
            {
                gd.addNumericField("Radius in Angstrom", radius_angstrom, 2);   
            }
            else
            {
                gd.addNumericField("Radius in Pixel", radius, 2);
            }
            gd.addCheckbox("accept radius", false);
            OvalRoi oroi=new OvalRoi(points[0].x-(int)radius/2,points[0].y-(int)radius/2,radius,radius);
            imp.setRoi(oroi);
            gd.showDialog();
            if(scale!=-1)
            {
                radius_angstrom=gd.getNextNumber();
                radius=radius_angstrom*imheight/(fov*scale);
                
            }
            else
            {
                radius=gd.getNextNumber();
            }
            baccepted=gd.getNextBoolean();
            
            if(gd.wasCanceled())
            {
                return false;
            }
        }
        return true;
        
    }
    
    double getScale(String unit)
    {
        if (unit.equals("nm"))
            return 10;
        else if (unit.equals("A"))
            return 1;
        else if (unit.equals("pm"))
            return 0.01;
        else if (unit.equals("um") ||  (unit.equals("µm")))
            return 10000;
        else
            return -1;           
    }
    
    OvalRoi optimizeRoi(int x, int y)
    {
        OvalRoi oroi = new OvalRoi(x,y,radius,radius);
        imp.setRoi(oroi);
        int bestx=x;
        int besty=y;
        double bestmean=imp.getStatistics().mean;
        int tolradius=3;
        //System.out.println("mean before: "+bestmean);
        for (int xs =-tolradius;xs<=tolradius;xs++)
        {
            for (int ys =-tolradius;ys<=tolradius;ys++)
            {
                oroi.setLocation(x+xs, y+ys);
                //imp.setRoi(oroi);
                ImageStatistics is = imp.getStatistics();
                double newmean=is.mean;
                if(newmean>bestmean)
                {
                    bestmean=newmean;
                    bestx=x+xs;
                    besty=y+ys;    
                }
                /*try {
                    TimeUnit.SECONDS.sleep(1);
                } catch (InterruptedException ex) {
                    Logger.getLogger(Calculate_Intensities.class.getName()).log(Level.SEVERE, null, ex);
                }*/
        
            }
        }
        //System.out.println("mean new: "+bestmean);
        if(!meanarray.add(bestmean))
            System.out.println("ERROR in adding mean");
        oroi.setLocation(bestx, besty);
        roim.addRoi(oroi);  
        return oroi;
        
    }
    
    double calculateMean(ArrayList<Double> al)
    {
        double ret=0;
        if(al.size()==0)
            return 0;
        
        for(Double d : al)
        {
            ret+=d;
        }
        ret/=al.size();
        return ret;   
    }
    
    double calculateStandardDeviation(ArrayList<Double> al)
    {
        double ret=0;
        if(al.size()==0)
            return 0;
        double mean=calculateMean(al);
        for(Double d : al)
        {
            ret=ret+(d-mean)*(d-mean);
        }
        ret/=al.size();
        ret=Math.sqrt(ret);
        return ret;   
    }

}
