package com.smallnn.input;

import static org.apache.commons.lang.SystemUtils.getUserHome;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FileUtil {

    public static List<File> find(File file) {
        return find(file, null);
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

    public static File tildeExpand(String path) {
        if (path.startsWith("~")) {
            path = path.replaceFirst("~", getUserHome().getAbsolutePath());
        }
        return new File(path);
    }
}
