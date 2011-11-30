package org.phymod.tools;

import java.util.Random;

public class Rand {
	/**
	 * creates an array of random integers
	 * 
	 * @param count number of random values to create
	 * @param max maximum value
	 * @param distinct if true each int will occur only once
	 * @return
	 */
	public static int[] createArray(int count, int max, boolean distinct) {
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
	
	/**
	 * return new integer random 
	 * 
	 * @param max value
	 * @return
	 */
	public static int createNew(int max) {
		Random rand = new Random();
		
		return rand.nextInt(max);
	}
}
