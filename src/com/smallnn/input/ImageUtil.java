package com.smallnn.input;

import static com.smallnn.AlgebraUtil.featureNormalize;
import static com.smallnn.input.FileUtil.filenames;

import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.vecmath.GMatrix;

import org.junit.Test;

import com.smallnn.AlgebraUtil.Normalized;

public class ImageUtil {
    public static final int SMALL_HEIGHT = 5;
    public static final int SMALL_WIDTH = 6;
    public static final int HEIGHT = 60;
    public static final int WIDTH = 72;

    public static final int IMAGE_SIZE = WIDTH * HEIGHT;
    
    public static class Resolution {
        public final int width;
        public final int height;
        public Resolution(int width, int height){
            this.width = width;
            this.height = height;
        }
    }

    public static long[] imageToLongArray(BufferedImage img){
        byte[] tmp = ((DataBufferByte) (img.getRaster().getDataBuffer())).getData();
        long[] result = new long[tmp.length];
        for(int i=0;i<tmp.length;i++){
            result[i] = tmp[i] & 0xffl;
        }
        return result;
    }
    
    public static long[] resizedImageToLongArray(BufferedImage resized) throws Exception {
        int[] tmp = ((DataBufferInt) (resized.getRaster().getDataBuffer())).getData();
        long[] result = new long[tmp.length];
        for(int i=0;i<tmp.length;i++){
            result[i] = tmp[i] & 0xffffffffl;
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
    
    public static long[] readResizeImage(BufferedImage img, Resolution orig, Resolution fin) throws Exception {
        if (fin.width == orig.width && fin.height == orig.height)
            return imageToLongArray(img);
        
        return resizedImageToLongArray(resize(img, fin.width, fin.height));
    }
    
    public static void printImageVals(BufferedImage src){
        BufferedImage img2 = resize(src, SMALL_WIDTH, SMALL_HEIGHT);
        long[] second = imageToLongArray(img2);
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
        BufferedImage resizedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        Graphics2D g = resizedImage.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

        g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.drawImage(image, 0, 0, width, height, null);
        g.dispose();
        return resizedImage;
    }

    int idx = 0;

    @Test
    public void testImageDisplay() throws Exception {
        JFrame frame = new JFrame("ImageUtil.testImageDisplay");
        frame.getContentPane().paint(resize(ImageIO.read(new File(filenames[1])), SMALL_WIDTH, SMALL_HEIGHT).getGraphics());
        frame.setSize(200, 200);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
}
