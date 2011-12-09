package com.smallnn.input;

import static com.smallnn.AlgebraUtil.featureNormalize;
import static com.smallnn.input.FileUtil.filenames;
import static com.smallnn.input.ImageUtil.HEIGHT;
import static com.smallnn.input.ImageUtil.IMAGE_SIZE;
import static com.smallnn.input.ImageUtil.SMALL_HEIGHT;
import static com.smallnn.input.ImageUtil.SMALL_WIDTH;
import static com.smallnn.input.ImageUtil.WIDTH;
import static com.smallnn.input.ImageUtil.imageToLongArray;
import static com.smallnn.input.ImageUtil.resize;

import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.vecmath.GMatrix;

import org.junit.Test;

import com.smallnn.AlgebraUtil.Normalized;

public class ImageUtilTest {

    @Test
    public void testGrayscale() throws Exception{
        BufferedImage img2 = resize(ImageIO.read(new File(filenames[2])), SMALL_WIDTH, SMALL_HEIGHT);
        BufferedImage img3 = resize(ImageIO.read(new File(filenames[3])), SMALL_WIDTH, SMALL_HEIGHT);
        long[] second = imageToLongArray(img2);
        long[] third = imageToLongArray(img3);
        for (long b: second){
            System.out.print(b+" ");
        }
        System.out.println();
        for (long b: third){
            System.out.print(b+" ");
        }
        System.out.println();
        
        double[] dSedondAndThird = new double[2*second.length+third.length];
        for(int i=0;i<second.length;i++){
            dSedondAndThird[i] = (double) second[i];
        }
        for(int i=0;i<second.length;i++){
            dSedondAndThird[second.length+i] = (double) second[i];
        }
        
        for(int i=0;i<third.length;i++){
            dSedondAndThird[2*second.length+i] = (double) third[i];
        }
        
        GMatrix mtx = new GMatrix(3, second.length, dSedondAndThird);
        Normalized nrm = featureNormalize(mtx);
        System.out.print(nrm.values);
    }
    
    
    @Test
    public void smallImageInit() throws Exception {
        BufferedImage img = ImageIO.read(new File(filenames[0]));
        Graphics smallImg = img.getGraphics().create(0, 0, 6, 5);

        JFrame frame = new JFrame("Display image");
        frame.getContentPane().print(smallImg);
        frame.setSize(10, 10);
        frame.setVisible(true);
    }
    
    @Test
    public void test() throws Exception {
        double[][] x = new double[filenames.length][IMAGE_SIZE];
        
        for (int l = 0; l < filenames.length; l++) {

            BufferedImage img = ImageIO.read(new File(filenames[l]));

            int height = img.getHeight();
            assert height == HEIGHT;
            for (int i = 0; i < height; i++) {
                int width = img.getWidth();
                assert width == WIDTH;
                for (int j = 0; j < width; j++) {
                    int pixelInx = i * height + j;
                    int rgb = img.getRGB(j, i);
                    x[l][pixelInx] = rgb;
                }
            }
        }
    }

}
