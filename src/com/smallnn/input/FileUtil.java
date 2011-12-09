package com.smallnn.input;

import static org.apache.commons.lang.SystemUtils.getUserHome;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingDeque;

import org.junit.Test;

public class FileUtil {

    public static List<File> find(File file) {
        return find(file, null);
    }
    
    @Test
    public void testFileList() throws Exception {
        
        for(File f : find(tildeExpand("~/References/barsDetect/training"), pngFilter)){
            System.out.println(f.getAbsolutePath());
        }
    }
    
    public static final FileFilter pngFilter  = new FileFilter() {

        @Override
        public boolean accept(File f) {
            return f.isFile() && f.getName().toLowerCase().endsWith("png");
        }
        
    };
    
    public static final FileFilter dirFilter  = new FileFilter() {

        @Override
        public boolean accept(File f) {
            return f.isDirectory();
        }
        
    };
    
    public static List<File> find(File file, FileFilter filter) {
        List<File> result = new ArrayList<File>();
        LinkedList<File> stack = new LinkedList<File>();
        stack.push(file);
        while (!stack.isEmpty()) {
            File f = stack.pop();
            if (filter == null || filter.accept(f)) {
                result.add(f);
            }

            if (f.isDirectory() && f.exists()) {
                stack.addAll(Arrays.asList(f.listFiles()));
            }
        }
        return result;
    }

    public final static File EOF = new File("EOF");

    public static BlockingQueue<File> findInThread(final File file, final FileFilter filter) {
        final BlockingQueue<File> result = new LinkedBlockingDeque<File>();
        Thread t = new Thread() {
            public void run() {
                Deque<File> stack = new LinkedList<File>();
                stack.push(file);
                while (!stack.isEmpty()) {
                    File f = stack.pop();
                    if (filter == null || filter.accept(f)) {
                        result.add(f);
                    }

                    if (f.isDirectory() && f.exists()) {
                        List<File> asList = new ArrayList<File>(Arrays.asList(f.listFiles()));
                        Collections.reverse(asList);
                        for (File file2 : asList) {
                            stack.addFirst(file2);
                        }
                    }
                }
                result.add(EOF);
            }
        };
        t.setDaemon(true);
        t.start();
        return result;
    }

    public static File tildeExpand(String path) {
        if (path.startsWith("~")) {
            path = path.replaceFirst("~", getUserHome().getAbsolutePath());
        }
        return new File(path);
    }

    static final String[] filenames = {

            "/Users/arsen/References/barsDetect/training/bars/72x60/61.scaled.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-041.2.72x60.png",
            
            "/Users/arsen/References/barsDetect/training/bars/72x60/replacementkillers_2005_ec_CC_16x9_240_english_2136_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-041.2.72x60.png",

            "/Users/arsen/References/barsDetect/training/bars/72x60/jeopardy_5301_test_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-033.72x60.png",

            "/Users/arsen/References/barsDetect/training/bars/72x60/harttohart_sd_63_las_no7241_test_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-034.72x60.png",

            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-001.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-060.72x60.png",

            "/Users/arsen/References/barsDetect/training/bars/72x60/61.40x40.to.72x60.png",

            "/Users/arsen/References/barsDetect/training/bars/72x60/61.ffmpeg.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/61.rescaled.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/badteacher_2011_rated_hd_16x9_185_2398_english_6207_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/burlesque_2010_hd_16x9_240_2398_english_6146_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/fogthe_unrated_2005_hd_16x9_235_2398_english_74_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/groove_2000_hd_16x9_178_2398_english_4814_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/socialnetworkthe_2010_hd_16x9_240_2398_english_5135_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/steve_dont_get_nun_1361_tv_2727.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/swat_1975_16_sd_4x3_133_25_english_1180_JPEG2000.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-002.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-003.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-004.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-005.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-006.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-007.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-008.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-009.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-010.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-011.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-012.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-013.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-014.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-015.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-016.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-017.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-018.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-019.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-020.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-021.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-022.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-023.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-024.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-025.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-026.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-027.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-028.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-029.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-030.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-031.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-032.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-035.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-036.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-037.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-038.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-039.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-040.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-041.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-042.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-043.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-044.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-045.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-046.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-047.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-048.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-049.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-050.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-051.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-052.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-053.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-054.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-055.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-056.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-057.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-058.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-059.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-061.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-062.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-063.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-064.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-065.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-066.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-067.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-068.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-069.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-070.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-071.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-072.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-073.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-074.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-075.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-076.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-077.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-078.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-079.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-080.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-081.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-082.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-083.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-084.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-085.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-088.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-089.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-086.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-087.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-090.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-091.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-092.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-093.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-094.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-095.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-096.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/tvproxy-bar-100.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-042.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-043.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-044.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-045.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-046.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-047.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-040.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-048.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-049.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-050.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-002.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-003.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-004.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-005.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-006.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-001.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-007.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-008.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-009.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-010.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-097.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-011.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-099.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-012.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-013.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/bars/72x60/tvproxy-bar-098.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-014.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-015.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-016.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-017.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-018.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-019.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-020.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-021.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-022.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-023.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-024.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-025.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-026.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-027.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-028.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-029.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-030.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-031.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-032.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-033.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-034.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-035.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-036.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-037.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-038.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-039.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-051.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-052.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-053.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-054.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-055.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-056.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-057.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-058.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-059.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-060.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-061.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-062.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-063.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-064.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-065.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-066.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-067.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-068.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-069.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-070.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-071.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-072.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-073.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-074.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-075.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-076.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-077.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-078.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-079.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-080.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-081.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-082.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-083.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-084.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-085.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-086.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-087.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-088.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-089.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-090.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-091.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-092.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-093.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-094.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-095.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-096.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-097.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-098.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-099.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/bar-100.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-001.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-001.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-002.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-002.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-003.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-003.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-004.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-004.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-005.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-005.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-006.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-006.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-007.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-007.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-008.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-008.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-009.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-009.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-010.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-010.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-011.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-011.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-012.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-012.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-013.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-013.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-014.2.72x60.png",
            "/Users/arsen/References/barsDetect/training/non-bars/72x60/non-bar-014.72x60.png" };

    public static final String[] testcases = { "/Users/arsen/References/barsDetect/test/bars/72x60/bar-001.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-100.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-002.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-003.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-004.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-005.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-006.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-007.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-008.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-009.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-010.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-011.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-012.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-013.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-014.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-015.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-016.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-017.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-018.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-019.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-020.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-021.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-022.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-023.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-024.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-025.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-026.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-027.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-028.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-029.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-030.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-031.72x60.png",
            "/Users/arsen/References/barsDetect/test/bars/72x60/bar-032.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-001.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-002.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-003.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-004.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-005.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-006.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-007.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-008.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-009.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-010.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-011.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-012.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-013.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-014.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-015.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-016.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-017.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-018.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-019.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-020.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-021.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-022.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-023.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-024.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-025.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-026.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-027.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-028.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-029.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-030.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-031.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-032.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-033.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-034.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-035.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-036.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-037.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-038.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-039.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-040.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-041.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-042.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-043.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-044.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-045.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-046.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-047.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-048.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-049.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-050.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-051.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-052.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-053.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-054.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-055.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-056.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-057.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-058.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-059.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-060.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-061.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-062.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-063.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-064.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-065.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-066.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-067.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-068.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-069.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-070.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-071.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-072.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-073.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-074.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-075.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-076.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-077.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-078.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-079.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-080.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-081.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-082.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-083.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-084.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-085.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-086.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-087.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-088.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-089.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-090.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-091.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-092.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-093.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-094.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-095.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-096.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-097.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-098.72x60.png",
            "/Users/arsen/References/barsDetect/test/non-bars/72x60/bar-099.72x60.png", };

}
