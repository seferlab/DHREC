//History Reconstruction Plot Generator
import java.util.*;
import java.io.File;
import java.io.IOException;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;


public class AxisInfo(){
    ArrayList xticks;
    ArrayList yticks;
    String xlabel;
    String ylabel; 
    String title;    

    public void AxisInfo(ArrayList<xtype> xticks, ArrayList<ytype> yticks, String xlabel, String ylabel, String title){
	 this.xticks = new ArrayList<xtype>(xticks);
         this.yticks = new ArrayList<ytype>(yticks);
         this.xlabel = xlabel;
         this.ylabel = ylabel;
         this.title = title;
    }
}

public class PlotGenerator(){
    private final String plotPath = "plots";
    private final String plotType = "fracs";
    private final String resultPref = "result";
    private final String tracePref = "traces";
    private final String realData = "real";
    private final float noise = 0.0;   
    private final int srate = 0;
    private final int startCount = 5;
    private final float timeFrac = 0.26;
    private final int timeCount = 3;    
    private final String evol = "static";  
    private final String smodel = "si";    
    private final String prob = "dis"; // "cont"   
    private final String inferMode = "Spreader"; //History     
    private final ArrayList<String> useStates = new ArrayList("i");    
    private ArrayList<String> curScores;    
    
    
    public static ArrayList<Float> getFracs(){
        return new ArrayList<Float>(Arrays.asList(0.1,0.2,0.3,0.5,0.75));
    }
    public static ArrayList<Float> getNoises(){
        return new ArrayList<Float>(Arrays.asList(0.1,0.2,0.3,0.5,0.7,0.9));
    }
    public static ArrayList<Integer> getStartcounts(){
        return new ArrayList<Integer>(Arrays.asList(1,3,5,7,9,10));
    }
    public static ArrayList<Integer> getTimecounts(){
        return new ArrayList<Integer>(Arrays.asList(1,2,3,4,5,10));   
    }


    /**
    * Generates heatmap given data and info and saves it as png
    * @param  x
    * @param  y
    * @param info
    * @param filename 
    * @return void
    * @see 
    */
    public static void genHeatmap(Map<String,ArrayList<Float>> x, Map<String,ArrayList<Float>> y, AxisInfo info,  String filename){
        StringBuilder tabsPath = new StringBuilder();
        tabsPath.append(this.plotPath);
        tabsPath.append("/");
        tabsPath.append(filename);
        String absPath = tabsPath.toString();  

        double[][] data = new double[][]{{3,2,3,4,5,6},
                                 {2,3,4,5,6,7},
                                 {3,4,5,6,7,6},
                                 {4,5,6,7,6,5}};
        HeatChart map = new HeatChart(data);
        map.setTitle("This is my heat chart title");
        map.setXAxisLabel("X Axis");
        map.setYAxisLabel("Y Axis");
        map.saveToFile(new File(absPath));   
    }

    
    /**
    * Generates scatter plot given x, y and info and saves it as png
    * @param  x
    * @param  y
    * @param info
    * @param filename 
    * @return void
    * @see 
    */
    public static void genPlot(Map<String,ArrayList<Float>> x, Map<String,ArrayList<Float>> y, AxisInfo info,  String filename){
        StringBuilder tabsPath = new StringBuilder();
        tabsPath.append(this.plotPath);
        tabsPath.append("/");
        tabsPath.append(filename);
        String absPath = tabsPath.toString();   
        
        dataset = new XYSeriesCollection();
        for(algo in x.keySet()){
            XYSeries data = new XYSeries(algo);
            for(int i=0; i < x.get(algo).size(); i++){
	        data.add(x.get(algo)[i], y.get(algo)[i])
            }
            dataset.addSeries(data);
	}
        
        final JFreeChart chart = ChartFactory.createScatterPlot(info.title, info.x, info.y, y, PlotOrientation.VERTICAL, true, true, false);
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        plot.setRenderer(renderer);
        
        File filename_png = new File(absPath);
        try {
           ChartUtilities.saveChartAsPNG(filename_png, chart, 980, 550);
        } catch (IOException ex) {
           throw new RuntimeException("Error saving a file",ex);
        }
                   
    }

    public static void main(String args[]){
        if (this.inferMode.equals("Spreader")){
	   this.curScores = ArrayList();
        } else if(this.inferMode.equals("Spreader")){
           this.curScores = ArrayList("KenTauCorrelation");
        }     
             
    }     

}
