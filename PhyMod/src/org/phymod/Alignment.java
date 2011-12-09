package org.phymod;

import java.util.LinkedList;

import org.phymod.io.FileLoader;
import org.phymod.tools.Rand;

public class Alignment {
	protected LinkedList<Taxa> taxa;
	public int nrPartitions;
	public int nrTaxa;
	
	public Alignment() { }
	
	public Alignment(Taxa[] s) {
		taxa = new LinkedList<Taxa>();
		for (int i = 0 ; i < s.length ; i++)
			taxa.addLast(s[i]);
		
		update();
	}
	
	public Alignment(LinkedList<Taxa> s) {
		taxa = s;
		update();
	}
	
	/**
	 * updates number of taxa and partitions per taxa
	 * @return this instance
	 */
	private Alignment update() {
		nrPartitions = taxa.peekFirst().partitions.length;
		nrTaxa = taxa.size();
		return this;
	}
	
	/**
	 * create alignment and partition files
	 * @param name
	 */
	public void toFile(String name) {
		FileLoader phylip = new FileLoader(name + ".phy");
		FileLoader part = new FileLoader(name + ".part");
		phylip.write(this.phylipToString());
		part.write(this.partToString());
	}	
	
	public String phylipToString() {
		String s = " " + taxa.size() + " " + taxa.get(0).length + "\n";
		for(int i = 0; i < taxa.size(); i++) {
			s += taxa.get(i).name + " ";
			for(int j = 0; j < taxa.get(0).partitions.length; j++) {
				s += taxa.get(i).partitions[j].data;
			}
			s += "\n";
		}
		return s;
	}
	
	public String partToString() {
		String s = "";
		int start = 1;
		
		for(int j = 0; j < taxa.get(0).partitions.length; j++) {
			int end = start + taxa.get(0).partitions[j].data.length() - 1;
			s += taxa.get(0).partitions[j].model + ", " + taxa.get(0).partitions[j].name + " = " + start + "-" + end + "\n";
			start = end + 1;
		}
		return s;
	}

	public String toString() {
		String s = "";
		for(int i = 0; i < taxa.size(); i++) {
			s += taxa.get(i).name + "\n";
			for(int j = 0; j < taxa.get(i).partitions.length; j++) {
				s += taxa.get(i).partitions[j].data + "\n";
			}
		}
		return s;
	}
	
	/**
	 * removes taxons with indices not given in t
	 * @param t list of taxons to keep
	 * @return
	 */
	public Alignment reduceToTaxa(int[] t) {
		System.out.println("Reducing number of taxa to " + t.length);
		LinkedList<Taxa> tmp = new LinkedList<Taxa>();
		for(int i = 0; i < t.length; i++)
			tmp.addLast(taxa.get(t[i]));
		this.taxa = tmp;
		
		return this.update();
	}
	
	/**
	 * remove partitions not given in p
	 * @param p
	 * @return
	 */
	public Alignment reduceToPartitions(int[] p) {
		System.out.println("Reducing number of partitions to " + p.length);
		LinkedList<Taxa> tmp = new LinkedList<Taxa>();
		
		for(int j = 0; j < taxa.size(); j++) {
			Partition[] parts = new Partition[p.length];
			for(int i = 0; i < p.length; i++)
				parts[i] = taxa.get(j).partitions[p[i]];
			
			tmp.addLast(new Taxa(taxa.get(j).name, parts));
		}
		this.taxa = tmp;
		return this.update();
	}
	
	/**
	 * randomly choose count partitions among all available ones
	 * 
	 * @param count number of partitions
	 * @param del if true distinct partitions are chosen
	 * @return
	 */
	public Alignment extractRand(int count, boolean del) {
		boolean redo;
		int[] r = new int[count];
		
		do {
			redo = false;
			r = Rand.createArray(count, this.nrPartitions, del);
		
			for(int j = 0; j < this.taxa.size(); j++) {
				String str = "";
				
				for(int i = 0; i < r.length; i++)
					str += this.taxa.get(j).partitions[r[i]].data;
				
				if(str.replaceAll("-", "").equals("") || str.replaceAll("X", "").equals("")) {
					redo = true;
					break;
				}
			}
		} while (del && redo);
		
		return this.reduceToPartitions(r);
	}
}
