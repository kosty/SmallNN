import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Panel;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JSplitPane;

import org.junit.Test;

public class ImageTest extends Panel {
    private static final int SMALL_HEIGHT = 5;
    private static final int SMALL_WIDTH = 6;
    private static final int HEIGHT = 60;
    private static final int WIDTH = 72;

    public static final int imageSize = WIDTH * HEIGHT;

    static final String[] filenames = { "/Users/arsen/References/misc/bars/61.cropped.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-041.2.72x60.png",
    /*
     * "/Users/arsen/References/misc/bars/72x60/61.40x40.to.72x60.png",
     * "/Users/arsen/References/misc/bars/72x60/61.rescaled.72x60.png",
     * "/Users/arsen/References/misc/bars/72x60/61.scaled.72x60.png",
     * "/Users/arsen/References/misc/bars/72x60/badteacher_2011_rated_hd_16x9_185_2398_english_6207_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/burlesque_2010_hd_16x9_240_2398_english_6146_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/fogthe_unrated_2005_hd_16x9_235_2398_english_74_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/groove_2000_hd_16x9_178_2398_english_4814_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/harttohart_sd_63_las_no7241_test_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/jeopardy_5301_test_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/replacementkillers_2005_ec_CC_16x9_240_english_2136_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/socialnetworkthe_2010_hd_16x9_240_2398_english_5135_JPEG2000.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/steve_dont_get_nun_1361_tv_2727.72x60.png"
     * ,
     * "/Users/arsen/References/misc/bars/72x60/swat_1975_16_sd_4x3_133_25_english_1180_JPEG2000.72x60.png"
     */
    };
    double[][] x = new double[filenames.length][imageSize];

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
        // for (int l = 0; l < filenames.length; l++) {
        // StringBuilder sb = new StringBuilder();
        // int cnt = 0;

        BufferedImage img = ImageIO.read(new File(filenames[0]));
        Graphics smallImg = img.getGraphics().create(0, 0, 6, 5);

        JFrame frame = new JFrame("Display image");
        frame.getContentPane().print(smallImg);
        frame.setSize(10, 10);
        frame.setVisible(true);

        // }
    }

    /*
     * public class ShowImage extends Panel { BufferedImage image;
     * 
     * public ShowImage() throws IOException {
     * System.out.println("Enter image name\n"); BufferedReader bf = new
     * BufferedReader(new InputStreamReader( System.in)); String imageName =
     * bf.readLine(); File input = new File(imageName); image =
     * ImageIO.read(input); }
     * 
     * public void paint(Graphics g) { g.drawImage(image, 0, 0, null); }
     * 
     * }
     */

    private BufferedImage image;

    private static BufferedImage resize(BufferedImage image, int width, int height) {
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

    public ImageTest(String filename) {
        try {
            image = resize(ImageIO.read(new File(filename)), SMALL_WIDTH, SMALL_HEIGHT);
        } catch (IOException ie) {
            ie.printStackTrace();
        }
    }

    private JSplitPane splitPane;

    public void paint(Graphics g) {
        g.drawImage(image, 0, 0, null);
    }

    static public void main(String args[]) throws Exception {
        JFrame frame = new JFrame("ShowImage.java");
        frame.getContentPane().add(new ImageTest(filenames[1]));
        frame.setSize(200, 200);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
}
