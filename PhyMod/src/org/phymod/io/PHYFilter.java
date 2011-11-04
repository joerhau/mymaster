package org.phymod.io;

import java.io.File;
import java.io.FilenameFilter;

public class PHYFilter implements FilenameFilter {
	
	@Override
	public boolean accept(File dir, String name) {
		if(name.indexOf("phy") > 0) return true;
		return false;
	}

}
