package com;

import static com.FileUtil.find;
import static com.FileUtil.pngFilter;
import static com.FileUtil.tildeExpand;
import static com.ImageUtil.SMALL_HEIGHT;
import static com.ImageUtil.SMALL_WIDTH;
import static com.ImageUtil.imageToGrayscale;
import static com.ImageUtil.resize;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.vecmath.GMatrix;

public class TrainDataUtil {
    
    public static final String dataPrefix = "~/References/barsDetect";
    public static final String TRAINING_DIR = "training";
    public static final String VALIDATION_DIR = "verify";
    public static final String TEST_DIR = "test";
    
    public static File getDataDir(String set, String cluster, String resolution){
        return tildeExpand(dataPrefix+"/"+set+"/"+cluster+"/"+resolution);
    }

    public static File getDataDir(String set){
        return tildeExpand(dataPrefix+"/"+set);
    }
    
    public static double[] longToDoubleArray(long[] arr){
        double[] result = new double[arr.length];
        for(int i=0;i<arr.length;i++){
            result[i] = (double)arr[i];
        }
        return result;
    }
    
    public static Data getMixedData() throws Exception{
        List<File> trainingFileSet = interleaveEntries();
        int inputLayerSize = SMALL_WIDTH*SMALL_HEIGHT;
        int classes = 2;
        int m=trainingFileSet.size()*40;
        double[] x = new double[m*inputLayerSize];
        double[] y = new double[m*classes];
        
        for(int i=0;i<m;i++){
            File trainFile = trainingFileSet.get(i%trainingFileSet.size());
            double[] tmparr = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(trainFile), SMALL_WIDTH, SMALL_HEIGHT))); 
            for (int j=0;j<inputLayerSize;j++){
                x[i*inputLayerSize+j] = tmparr[j];
            }
            y[i*classes]   = i%2 == 0 ? 1. : 0.;
            y[i*classes+1] = i%2 == 0 ? 0. : 1.;
        }
        
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }
    
    public static GMatrix readData(File f) throws Exception{
        double[] data = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(f), SMALL_WIDTH, SMALL_HEIGHT)));
        return new GMatrix(1, SMALL_WIDTH*SMALL_HEIGHT, data);
    }
    
    
    public static List<File> interleaveEntries(){
        List<File> bars = find(getDataDir(TRAINING_DIR, "bars", "72x60"), pngFilter);
        List<File> nonBars = find(getDataDir(TRAINING_DIR, "non-bars", "72x60"), pngFilter);
        List<File> interleaved = new ArrayList<File>();
        for(int i=0;i<bars.size();i++){
            interleaved.add(bars.get(i));
            interleaved.add(nonBars.get(i%nonBars.size()));
        }
        return interleaved;
    }
    
    public static List<File> listSetFiles(String set) throws Exception {
        return find(getDataDir(set), pngFilter);
    }
    
    public static class Data {
        
        public final GMatrix x;
        public final GMatrix y;

        public Data(GMatrix x, GMatrix y) {
            this.x = x;
            this.y = y;
        }
        
    }
    
}
