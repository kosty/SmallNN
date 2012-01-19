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
package com.smallnn.input;

import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class ImageUtil {
    
    public static class Resolution {
        public final int width;
        public final int height;
        public Resolution(int width, int height){
            this.width = width;
            this.height = height;
        }
    }

    public static double[] imageToDoubleArray(BufferedImage img){
        byte[] tmp = ((DataBufferByte) (img.getRaster().getDataBuffer())).getData();
        double[] result = new double[tmp.length];
        for(int i=0;i<tmp.length;i++){
            long l = tmp[i] & 0xffl;
            result[i] = (double) l;
        }
        return result;
    }
    
    public static double[] resizedImageToDoubleArray(BufferedImage resized) throws Exception {
        int[] tmp = ((DataBufferInt) (resized.getRaster().getDataBuffer())).getData();
        double[] result = new double[tmp.length];
        for(int i=0;i<tmp.length;i++){
            long l = tmp[i] & 0xffffffffl;
            result[i] = (double) l;
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
    
    public static BufferedImage zealousSubsCrop(BufferedImage img){
        int initialHeight = img.getHeight();
        int initialWidth  = img.getWidth();
        int newHeight = initialHeight/2;
        int newWidth  = initialWidth/2;
        return img.getSubimage(newWidth/2, newHeight, newWidth, newHeight);
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

}
