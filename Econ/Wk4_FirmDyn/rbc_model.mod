// RBC Model with additively separable preferences over labor
//
// Thomas Winberry, July 19th, 2019

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables
//----------------------------------------------------------------

var 

// Allocation variables 
y c k i l

// Prices
r w 

// Productivity shock
z

// Auxiliary variables for analysis
logY logC logI logL logW r_annual;

//----------------------------------------------------------------
// 2. Exogenous variables
//----------------------------------------------------------------

varexo e;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Utility function
bbeta eeis eeta cchi l_ss

// Technology
aalpha ddelta rrho ssigma;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences 
bbeta   	= 0.99;
eeis 		= 1;
l_ss 		= 1 / 3;
eeta 		= 1 / 2;

// Technology
aalpha 	 	= 1 / 3;
ddelta  	= 0.025;
rrho 		= 0.95;
ssigma  	= 0.007;

//----------------------------------------------------------------
// 5. Model
//----------------------------------------------------------------

model; 

	// 1. Euler equation for capital
	c ^ (-1 / eeis)				= bbeta * (c(+1) ^ (-1 / eeis)) * (1 - ddelta + r(+1));
	
	// 2. First order equation for labor-leisure
	w * (c ^ (-1 / eeis))			= cchi * (l ^ (1 / eeta));
	
	// 3. Production function
	y					= exp(z) * (k(-1) ^ aalpha) * (l ^ (1 - aalpha));
	
	// 4. Wage rate
	w					= (1 - aalpha) * y / l;
	
	// 5. Rental rate on capital
	r					= aalpha * y / k(-1);
	
	// 6. Evolution of capital stock
	k					= (1 - ddelta) * k(-1) + i;
	
	// 7. Output market clearing
	y					= c + i;
	
	// 8. Law of motion for productivity
	z					= rrho * z(-1) + e;
	
	// 9. Definition of logY
	logY					= 100 * log(y);
	
	// 10. Definition of logC
	logC					= 100 * log(c);
	
	// 11. Definition of logI
	logI					= 100 * log(i);
	
	// 12. Definition of logL
	logL					= 100 * log(l);
	
	// 13. Definition of logW
	logW					= 100 * log(w);
	
	// 14. Definition of r_annual
	r_annual				= 400 * r;
	
end;

//----------------------------------------------------------------
// 6. Steady State
//----------------------------------------------------------------

steady_state_model;

	r		= (1 / bbeta) - (1 - ddelta);
	l		= l_ss;
	k		= (aalpha * (l ^ (1 - aalpha)) / r) ^ (1 / (1 - aalpha));
	y		= (k ^ aalpha) * (l ^ (1 - aalpha));
	i		= ddelta * k;
	c		= y - i;
	z		= 0;
	w		= (1 - aalpha) * y / l;
	cchi		= w * (c ^ (-1 / eeis)) / (l ^ (1 / eeta));
	
	logY		= 100 * log(y);
	logC		= 100 * log(c);
	logI		= 100 * log(i);
	logL		= 100 * log(l);
	logW		= 100 * log(w);
	r_annual	= 400 * r;
	
end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------

shocks;
	var e = ssigma ^ 2;
end;

steady;

check;

stoch_simul(hp_filter = 1600, irf = 80, order = 1) logY logC logI logL logW r_annual;
