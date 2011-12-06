package com.smallnn;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class PlotTest extends JPanel {
    double[] data = null;
    final int PAD = 20;

    public PlotTest(double[] data) {
        this.data = data;
    }

    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        int w = getWidth();
        int h = getHeight();
        g2.drawLine(PAD, PAD, PAD, h - PAD);
        g2.drawLine(PAD, h - PAD, w - PAD, h - PAD);
        double xScale = (w - 2 * PAD) / (data.length + 1);
        double maxValue = 100.0;
        double yScale = (h - 2 * PAD) / maxValue;
        // The origin location.
        int x0 = PAD;
        int y0 = h - PAD;
        g2.setPaint(Color.red);
        for (int j = 0; j < data.length; j++) {
            int x = x0 + (int) (xScale * (j + 1));
            int y = y0 - (int) (yScale * data[j]);
            g2.fillOval(x - 2, y - 2, 4, 4);
        }
    }
    
    public static void main(String[] args) {
        double[] data = {16.07526661105712, 11.717426215964569, 8.741818183278063, 6.638388155235342, 5.125496480739991, 4.03330414801575, 3.2443085252871358, 2.674274362584106, 2.26242600446126, 1.9648646536953007}; // = { 25, 60, 42, 75 };
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(new PlotTest(data));
        f.setSize(400,400);
        f.setLocation(200,200);
        f.setVisible(true);
    }
    
    
}
