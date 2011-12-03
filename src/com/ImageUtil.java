package com;

import static com.FileUtil.filenames;
import static com.SingleLayerNetwork.featureNormalize;

import java.awt.Canvas;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.vecmath.GMatrix;

import org.junit.Test;

import com.SingleLayerNetwork.Normalized;

public class ImageUtil {
    public static final int SMALL_HEIGHT = 5;
    public static final int SMALL_WIDTH = 6;
    public static final int HEIGHT = 60;
    public static final int WIDTH = 72;

    public static final int IMAGE_SIZE = WIDTH * HEIGHT;
    
    double[][] x = new double[filenames.length][IMAGE_SIZE];

    @Test
    public void test() throws Exception {
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

    @Test
    public void smallImageInit() throws Exception {
        BufferedImage img = ImageIO.read(new File(filenames[0]));
        Graphics smallImg = img.getGraphics().create(0, 0, 6, 5);

        JFrame frame = new JFrame("Display image");
        frame.getContentPane().print(smallImg);
        frame.setSize(10, 10);
        frame.setVisible(true);
    }

    private BufferedImage image;
    
    public static long[] imageToGrayscale(BufferedImage img){
        int[] tmp = ((DataBufferInt) (img.getRaster().getDataBuffer())).getData();
        long[] result = new long[tmp.length];
        for(int i=0;i<tmp.length;i++){
            result[i] = tmp[i] & 0xffffffl;
        }
        return result;
    }
    
    public static BufferedImage readImage(String imageFile)  {
        try {
            return ImageIO.read(new File(imageFile));
        } catch (IOException e) {
            System.err.println(imageFile);
            return null;
        }
    }
    
    @Test
    public void testGrayscale() throws Exception{
//        printImageVals(ImageIO.read(new File(filenames[0])));
//        printImageVals(ImageIO.read(new File(filenames[1])));
//        printImageVals(ImageIO.read(new File(filenames[2])));
//        printImageVals(ImageIO.read(new File(filenames[3])));
        BufferedImage img2 = resize(ImageIO.read(new File(filenames[2])), SMALL_WIDTH, SMALL_HEIGHT);
        BufferedImage img3 = resize(ImageIO.read(new File(filenames[3])), SMALL_WIDTH, SMALL_HEIGHT);
        long[] second = imageToGrayscale(img2);
        long[] third = imageToGrayscale(img3);
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
    
    public static void printImageVals(BufferedImage src){
        BufferedImage img2 = resize(src, SMALL_WIDTH, SMALL_HEIGHT);
        long[] second = imageToGrayscale(img2);
        for (long b: second){
            System.out.print(b+" ");
        }
        System.out.println();
        double[] dsecond = new double[second.length];
        for(int i=0;i<second.length;i++){
            dsecond[i] = (double) second[i];
        }
        GMatrix mtx_second = new GMatrix(1, dsecond.length, dsecond);
        Normalized nrm_second = featureNormalize(mtx_second);
        System.out.print(nrm_second.values);
    }
    
    public static long[] imageToRGB(BufferedImage img){
        int imgSize = getImageSize(img);
        long[] result = new long[imgSize];
        
        int height = img.getHeight();
        assert height == 5;
        for (int i = 0; i < height; i++) {
            int width = img.getWidth();
            assert width == 6;
            for (int j = 0; j < width; j++) {
                int pixelInx = i * height + j;
                result[pixelInx] = (long) (img.getRGB(j, i) & 0xffffffl);
            }
        }
        return result;
    }

    private static int getImageSize(BufferedImage img) {
        return img.getWidth()*img.getHeight();
    }

    public static BufferedImage resize(BufferedImage image, int width, int height) {
        BufferedImage resizedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = resizedImage.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

        g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.drawImage(image, 0, 0, width, height, null);
        g.dispose();
        return resizedImage;
    }

    int idx = 0;

    public void paint(Graphics g) {
        g.drawImage(image, 0, 0, null);
    }

    @Test
    public void testImageDisplay() throws Exception {
        JFrame frame = new JFrame("ImageUtil.testImageDisplay");
        frame.getContentPane().paint(resize(ImageIO.read(new File(filenames[1])), SMALL_WIDTH, SMALL_HEIGHT).getGraphics());
        frame.setSize(200, 200);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
    
    public static void main(String... args) throws Exception{
        
        JFrame frame = new JFrame("Image Demo");
        frame.setBounds(50, 50, 100, 100);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        List<BufferedImage> imgs = new ArrayList<BufferedImage>();
        for (String imageFile : filenames){
            imgs.add(resize(ImageIO.read(new File(imageFile)), SMALL_WIDTH, SMALL_HEIGHT));
        }
        frame.getContentPane().add(new GCanvas(imgs));
        frame.setVisible(true);
    }

    static class GCanvas extends Canvas {
        
        final List<BufferedImage> images;
        
        public GCanvas(List<BufferedImage> filenames){
            this.images = filenames;
        }
        
        public void paint(Graphics g) {
            int xoff = 0;
            for (BufferedImage img : images) {
                g.drawImage(img, xoff, 0, this);
                xoff += img.getWidth();
            }
        }
    }
}
