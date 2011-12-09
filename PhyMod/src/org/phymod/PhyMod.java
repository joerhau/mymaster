package org.phymod;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.io.comparator.NameFileComparator;

import org.phymod.io.FileLoader;
import org.phymod.io.PHYFilter;
import org.phymod.io.PhylipLoader;
import org.phymod.io.PhylipPartitionLoader;
import org.phymod.tools.Rand;

public class PhyMod {
	public static void main(String[] args) {
		if(args.length == 0 || args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("--help")) {
			printHelp();
			System.exit(0);
		}
		
		parseArgs(args);
	}
	
	/**
	 * check given arguments
	 * 
	 * @param args
	 * @throws IOException
	 */
	private static void parseArgs(String[] args) {
		Command c ;
		if(args[0].equalsIgnoreCase("--glue") || args[0].equalsIgnoreCase("glue"))
			c = Command.GLUE;
		else if(args[0].equalsIgnoreCase("--assign") || args[0].equalsIgnoreCase("assign"))
			c = Command.ASSIGN;
		else if(args[0].equalsIgnoreCase("--extract") || args[0].equalsIgnoreCase("extract"))
			c=Command.EXTRACT;
		else
			c = Command.valueOf(args[0]);
		
		switch(c) {
		case EXTRACT: extract(args); break;
		
		case ASSIGN: {
			if(args.length < 2) {
				printHelp();
				System.exit(0);
			}
			assign(args);
		} break;
		
		case GLUE: {
			if(args.length < 2) {
				printHelp();
				System.exit(0);
			}
			File f = new File(args[1]);
			glue(f);
		} break;
		default: {
			System.out.println("Unknown Command passed");
			System.exit(0);
		}
		}
	}
	
	private static void extract(String[] args) {
		boolean del = false;
		int count = 0;
		boolean p = false;
		boolean m = true;
		boolean s = false;
		boolean r = false;
		int scount = 0;
		String tmp = "";
		String outfile = "";
		int[] reduce, remove;
		
		Alignment phy = new PhylipLoader(args[args.length - 1], args[args.length - 2]).getAlignment();
		String workdir = new File(new File(args[args.length - 1]).getAbsolutePath()).getParent();
		
		for(int i = 1; i < args.length; i++) {
			if(args[i].substring(0,1).equals("-")) {
//				choose random partitions
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
				else if(args[i].substring(1,2).equals("r")) {
					r = true;
					tmp = args[++i];
				}
//				name given
				else if(args[i].substring(1,2).equals("n")) outfile = args[++i];
//				do until no duplicates and non-determined species
				else if(args[i].substring(1,2).equals("d")) del = true;
			} 
		}
		
		// partitions to extract specified
		if(p) {
			String[] str = tmp.split(",");
			reduce = new int[str.length];
			for(int i = 0; i < str.length; i++)
				if(str[i].replaceAll(",", "") != "")
					reduce[i] = Integer.parseInt(str[i].replaceAll(",", ""));
			
			phy.reduceToPartitions(reduce);
		// remove the followong partitions
		} else if (r) {
			String[] str = tmp.split(",");
			remove = new int[str.length];
			for(int i = 0; i < str.length; i++)
				if(str[i].replaceAll(",", "") != "")
					remove[i] = Integer.parseInt(str[i].replaceAll(",", ""));
			
			reduce = new int[phy.nrPartitions - remove.length];
			
			for(int i = 0, j = 0; i < phy.nrPartitions; i++) {
				if(!intArrayContains(remove, i)) {
					reduce[j] = i;
					j++;
				}
			}
			
			phy.reduceToPartitions(reduce);
		// randomly extract some partitions
		} else if (m) {
			if(del && count > phy.nrPartitions) {
				System.out.println("Number of partitions to extract is larger than available ones. \"-d\" would result in an infinite Loop. Thus we exit here...");
				System.exit(0);
			}
			phy.extractRand(count, del);
		} 
		
		if(s) phy.reduceToTaxa(Rand.createArray(scount, phy.taxa.size(), true));
			
		if(!outfile.equals("")) phy.toFile(workdir + "/" + outfile);
	}
	
	private static boolean intArrayContains(int[] a, int c) {
		boolean ret = false;
		for(int i = 0; i < a.length; i++)
			if(a[i] == c) ret = true;
		return ret;
	}
	
	private static void assign(String[] args) {
		File assignmentFile = new File(args[1]);
		PhylipPartitionLoader part = new PhylipPartitionLoader(args[1]);
		
		FileLoader l = new FileLoader(assignmentFile);
		l.open();
		if(part.models.length == l.l.size()) {
			for(int i = 0; i < l.l.size(); i++) {
				part.models[i] = l.l.get(i);
			}
			System.out.println("going to create " + assignmentFile.getAbsolutePath() + ".part" + " using optimal model assignment");
			FileLoader write = new FileLoader(assignmentFile.getAbsolutePath() + ".part");
			write.write(part.toString());
		}
	}
	
	private static void glue(File f) {
		File[] files = f.listFiles(new PHYFilter());
		Arrays.sort(files, NameFileComparator.NAME_COMPARATOR);
		Alignment[] toGlue = new Alignment[files.length];
		
		int lastDot = files[0].getAbsolutePath().lastIndexOf(".");
		toGlue[0] = new PhylipLoader(files[0].getAbsolutePath(), files[0].getAbsolutePath().substring(0, lastDot) + ".part").getAlignment();
		LinkedList<Taxa> spec = new LinkedList<Taxa>();
		
		for(int j = 0; j < toGlue[0].taxa.size(); j++) 
			spec.add(j, toGlue[0].taxa.get(j));
				
		for(int i = 1; i < files.length; i++) {
			if(!files[i].getName().contains("synthetic")){
		    lastDot = files[i].getAbsolutePath().lastIndexOf(".");
			toGlue[i] = new PhylipLoader(files[i].getAbsolutePath(), files[i].getAbsolutePath().substring(0, lastDot) + ".part").getAlignment();
			
				for(int j = 0; j < toGlue[i].taxa.size(); j++) 
					spec.get(j).addPartition(toGlue[i].taxa.get(j).partitions[0]);
			}
		}
		
		Alignment o = new Alignment(spec);
		System.out.println("writing files " + f.getAbsolutePath() + "/synthetic.*");
		o.toFile(f.getAbsolutePath() + "/synthetic");
	}
	
	/**
	 * print usage
	 * 
	 */
	private static void printHelp() {
		System.out.println("TODO Usage: PhyMod COMMAND [Options] [Arguments]\n");
		System.out.println("\tEXTRACT <Partitionfile> <Phylipfile>");
		System.out.println("\textracting specified part of the alignment");
		System.out.println("\t-m\tRandomly extract some of the partitions included in PHYLIPFILE");
		System.out.println("\t\t\"-m count\":\tNumber of Partitions to extract");
		System.out.println("\t\t\t\tif (count > #partitions in PHYLIPFILE) some partitions may occure more than once)\n");
		System.out.println("\t-d\tFind combination, so that (1) there are no duplicates and (2) there are no");
		System.out.println("\t\tundetermined species.\n");
		System.out.println("\t-n\tSpecify names of outputfiles.");
		System.out.println("\t\t\"-n test\":\tOutput will be written to \"test.phy\" and \"test.part\".\n");
		System.out.println("\t-p\tSpecify partitions to extract.");
		System.out.println("\t\t\"-p 0,5,1\":\tPartitions 0, 5 and 1 will be included in output in exactly that order.\n");
		System.out.println("\t-s\tSpecify number of Species to extract. Species will be randomly chosen");
		System.out.println("\t\t\"-s scount\":\tNumber of Speciees.\n");
		System.out.println("\t-r\tSpecify Partition to remove");
		System.out.println("\t\t\"-r partition_number\":\tThis partition will be deleted from the input alignment.\n");
		
		System.out.println("\tGLUE <Folder>");
		System.out.println("\tConcatenating all phy file within the given directory into a new large file.\n");
		
		System.out.println("\tASSIGN");
		System.out.println("\tUse given substitution models to create a new partition file for the given dataset\n");
	}
}
