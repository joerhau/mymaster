package org.phymod.io;

import java.util.LinkedList;

import org.phymod.Alignment;
import org.phymod.Taxa;

public abstract class AlignmentLoader {
	protected int nrPartitions;
	protected int nrSpecies;
	protected LinkedList<Taxa> spec;
	
	public Alignment alignment;
	
	public AlignmentLoader() {
		spec = new LinkedList<Taxa>();
	};
	
	public abstract void init();
	
	public abstract Alignment getAlignment();
}
