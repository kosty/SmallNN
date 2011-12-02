package com;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;



public class ChartUtil {
    
    public static void showChart(double[] trainingCosts) throws IOException {
        XYSeries series = new XYSeries("Cost of predictions");
        for(int i=0;i<trainingCosts.length;i++){
            series.add(i, trainingCosts[i]);
        }
        //         Add the series to your data set
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        //         Generate the graph
        JFreeChart chart = ChartFactory.createXYLineChart("Neural network training cost", "iterations", "cost", dataset, PlotOrientation.VERTICAL, true, true, false);
        ChartFrame frame = new ChartFrame("First", chart);
        frame.pack();
        frame.setVisible(true);
    }
}
