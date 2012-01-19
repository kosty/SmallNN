/**
 *  This file is part of SmallNN, a small neural network implementation
 *  Copyright (C) 2011, 2012 Arsen Kostenko <arsen.kostenko@gmail.com>
 *     
 *  SmallNN is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SmallNN is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SmallNN.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.smallnn.output;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;



public class ChartUtil {
    
    public static void showChart(Double[] trainingCosts) throws IOException {
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
