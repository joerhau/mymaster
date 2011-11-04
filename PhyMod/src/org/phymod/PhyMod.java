package org.phymod;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.io.comparator.NameFileComparator;

import org.phymod.io.PHYFilter;
import org.phymod.io.PartitionLoader;

public class PhyMod {
	private static PhylipAlignment phy;
	private static PartitionLoader part;
	public static boolean del = false;

	public static void main(String[] args) throws IOException {
		int count = 0;
		boolean p = false;
		boolean m = true;
		boolean s = false;
		boolean glue = false;
		int scount = 0;
		String tmp = "";
		String outfile = "";
		PhylipAlignment out;
		int[] r;
		File glueFolder = null;
		PhylipAlignment[] toGlue;
		
		if(args.length < 2) { 
			printHelp();
			System.exit(1);
		}
		
		for(int i = 0; i < args.length; i++) {
			if(args[i].substring(0,2).equals("--")) {
				if(args[i].substring(2,6).equals("glue")) {
					glueFolder = new File(args[++i]);
					glue = true;
					break;
				}
			}
			else if(args[i].substring(0,1).equals("-")) {
//				chose random partitions
				if(args[i].substring(1,2).equals("m")) {
					count = Integer.parseInt(args[++i]);
					m = true;
				}
//				partitions specified
				else if(args[i].substring(1,2).equals("p")) {
					p = true;
					tmp = args[++i];
				}
//				number of species to extract
				else if(args[i].substring(1,2).equals("s")) {
					s = true;
					scount = Integer.parseInt(args[++i]);
				}
//				name given
				else if(args[i].substring(1,2).equals("n")) outfile = args[++i];
//				do until no duplicates and non-determined species
				else if(args[i].substring(1,2).equals("d")) del = true;
			}
			else {
				part = new PartitionLoader(args[i]);
				phy = new PhylipAlignment(args[i++], part);
			}
		}
		
		if(glue) {
			File[] files = glueFolder.listFiles(new PHYFilter());
			Arrays.sort(files, NameFileComparator.NAME_COMPARATOR);
			toGlue = new PhylipAlignment[files.length];
			int lastDot = files[0].getCanonicalPath().lastIndexOf(".");
			toGlue[0] = new PhylipAlignment(files[0].getCanonicalPath(), new PartitionLoader(files[0].getCanonicalPath().substring(0, lastDot) + ".part"));
			Species[] spec = new Species[toGlue[0].spec.length];
			
			for(int j = 0; j < toGlue[0].spec.length; j++) 
				spec[j] = toGlue[0].spec[j];
					
			for(int i = 1; i < files.length; i++) {
				if(!files[i].getName().contains("synthetic")){
			    lastDot = files[i].getCanonicalPath().lastIndexOf(".");
				toGlue[i] = new PhylipAlignment(files[i].getCanonicalPath(), new PartitionLoader(files[i].getCanonicalPath().substring(0, lastDot) + ".part"));
				
					for(int j = 0; j < toGlue[i].spec.length; j++) 
						spec[j].addPartition(toGlue[i].spec[j].partitions[0]);
				}
			}
			
			out = new PhylipAlignment(spec);
			System.out.println("writing files " + glueFolder.getCanonicalPath() + "/synthetic.*");
			out.toFile(glueFolder.getCanonicalPath() + "/synthetic");
		}
		else {
			if(p) {
				String[] str = tmp.split(",");
				r = new int[str.length];
				for(int i = 0; i < str.length; i++) {
					if(str[i].replaceAll(",", "") != "")
						r[i] = Integer.parseInt(str[i].replaceAll(",", ""));
				}
			} else if (m) {
				if(del && count > phy.nrPartitions) {
					System.out.println("Number of partitions to extract is larger than available ones. \"-d\" would result in an infinite Loop. Thus we exit here...");
					System.exit(0);
				}
			    r = createRandPart(count);
			} else {
				r = new int[phy.nrPartitions];
				for(int i = 0; i < r.length; i++)
					r[i] = i;
			}
			if(s) out = extract(r, createRandSpec(scount));
			else out = extract(r);
			
			if(!outfile.equals("")) {
				out.toFile(outfile);
			}
			else System.out.println(out.phylipToString());
		}
	}
	
	public static PhylipAlignment extract(int[] part) {
		int[] pass = new int[phy.spec.length];
		for(int i = 0; i < phy.spec.length; i++)
			pass[i] = i;
		return extract(part, pass);
	}
	
	
	public static PhylipAlignment extract(int[] part, int[] spec) {
		Species[] s = new Species[spec.length];
	    String s1 = "", s2 = "";
	    for(int i = 0; i < part.length; i++)
	    	s1 += part[i] + " ";
	    
	    for(int i = 0; i< spec.length; i++)
	    	s2 += spec[i] + " ";
	    
	    System.out.print("Extracting partitions " + s1);
	    System.out.print("from Species " + s2 + "\n");
	    
	    for (int i = 0; i < spec.length; i++) {
			s[i] = new Species(phy.spec[spec[i]].name, part.length);
		}
			
		for (int i = 0; i < spec.length; i++) {
			for(int j = 0; j < part.length; j++) {
				s[i].replace(phy.spec[spec[i]], part[j], j);
			}
		}
		
		return new PhylipAlignment(s);
	}
	
	
	private static int[] createRand(int count, int max, boolean distinct) {
		Random rand = new Random();
		int[] r = new int[count];
		boolean rep;
		
		for(int i = 0; i < count; i++) {
	    	int tmp;
	    	do {
	    		rep = false;
	    		tmp = rand.nextInt(max);
	    		for(int j = 0; j < i; j++) {
	    			if(r[j] == tmp) {
	    				rep = true;
	    				break;
	    			}
	    		}
	    	} while(distinct && rep);
	    	
	    	r[i] = tmp;
	    }
		return r;
	}
	
	private static int[] createRandSpec(int count) {
		return createRand(count, phy.spec.length, true);
	}
	
	private static int[] createRandPart(int count) {
		boolean redo;
		int[] r = new int[count];
		
		do {
			redo = false;
			
			r = createRand(count, phy.nrPartitions, del);
		
			for(int j = 0; j < phy.spec.length; j++) {
				String str = "";
				
				for(int i = 0; i < r.length; i++)
					str += phy.spec[j].partitions[r[i]].data;
				
				if(str.replaceAll("-", "").equals("")) {
					redo = true;
					break;
				}
			}
		} while (del && redo);
		
		return r;
	}
	
	private static void printHelp() {
		System.out.println("Usage: PhyMod [Options] PHYLIPFILE PARTITIONFILE\n");
		System.out.println("\t-m\tRandomly extract some of the partitions included in PHYLIPFILE\n");
		System.out.println("\t\t\"-m count\":\tNumber of Partitions to extract");
		System.out.println("\t\t\t\tif (count > #partitions in PHYLIPFILE) some partitions may occure more than once)\n");
		System.out.println("\t-d\tFind combination, so that (1) there are no duplicates and (2) there are no");
		System.out.println("\t\tundetermined species.\n");
		System.out.println("\t-n\tSpecify names of outputfiles.\n");
		System.out.println("\t\t\"-n test\":\tOutput will be written to \"test.phylip\" and \"test.part\".\n");
		System.out.println("\t-p\tSpecify partitions to extract.\n");
		System.out.println("\t\t\"-p 0,5,1\":\tPartitions 0, 5 and 1 will be included in output in exactly that order.\n");
		System.out.println("\t-s\tSpecify number of Species to extract. Species will be randomly chosen\n");
		System.out.println("\t\t\"-s scount\":\tNumber of Speciees.\n");
	}
}
