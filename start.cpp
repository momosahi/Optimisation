//----------------- declarations ---------------------

// labels
int hot   =1; /* fenetre 1 */
int cold  =2; /* fenetre 2 */
int radb  =3;
int radt  =4;
int wall  =5;
int flowin  =6;
int flowout  =7;

// geometry
real xmin = 0;
real xmax = 5;
real ymin = 0;
real ymax = 4;

real rbmin = 2.5;
real rbmax = 3.5;
real rtmin = 0.5;
real rtmax = 1.5;

real ftmin = 3;
real ftmax = 3.5;
real fbmin = 1;
real fbmax = 1.5;

real hotmin = 1;
real hotmax = 2;
real coldmin = 2;
real coldmax = 3;

real cxmin = 2;
real cxmax = 2.1;
real cymin = 2;
real cymax = 3.5;

real c2xmin = 2.;
real c2xmax = 4;
real c2ymin = 1.5;
real c2ymax = 1.6;

// mesh
int size = 5;


// physics
real temphot = 14;
real tempcold = 10;
real coeffdiff = 0.1;

real target = 19;

//----------------- definition points de controle flux -----------

real[int] paramOpt(3);

paramOpt[0] = 8; //8
paramOpt[1] = 12; //12
paramOpt[2] = 0.5; //0.5

//-------------- définition géométrie bords ---------------------

border bottom1(t=xmin,fbmin) { x=t; y=ymin ;label=wall;};
border bottom2(t=fbmin,fbmax) { x=t; y=ymin ;label=flowout;};
border bottom3(t=fbmax,rbmin) { x=t; y=ymin ;label=wall;};
border bottom4(t=rbmin,rbmax) { x=t; y=ymin ;label=radb;};
border bottom5(t=rbmax,xmax) { x=t; y=ymin ;label=wall;};

border rightb(t=ymin,coldmin) { x=xmax; y=t ;label=wall;};
border rightc(t=coldmin,coldmax) { x=xmax; y=t ;label=cold;};
border rightt(t=coldmax,ymax) { x=xmax; y=t ;label=wall;};

border top1(t=xmax,ftmax){ x=t; y=ymax; label=wall;};
border top2(t=ftmax,ftmin){ x=t; y=ymax; label=flowin;};
border top3(t=ftmin,rtmax){ x=t; y=ymax; label=wall;};
border top4(t=rtmax,rtmin){ x=t; y=ymax; label=radt;};
border top5(t=rtmin,xmin){ x=t; y=ymax; label=wall;};

border leftt(t=ymax,hotmax){ x=xmin; y=t ;label=wall;};
border lefth(t=hotmax,hotmin){ x=xmin; y=t ;label=hot;};
border leftb(t=hotmin,ymin){ x=xmin; y=t ;label=wall;};

border centerb(t=cxmax,cxmin){ x=t; y=cymin ;label=wall;};
border centerr(t=cymax,cymin){ x=cxmax; y=t ;label=wall;};
border centert(t=cxmin,cxmax){ x=t; y=cymax ;label=wall;};
border centerl(t=cymin,cymax){ x=cxmin; y=t ;label=wall;};

border center2b(t=c2xmax,c2xmin){ x=t; y=c2ymin ;label=wall;};
border center2r(t=c2ymax,c2ymin){ x=c2xmax; y=t ;label=wall;};
border center2t(t=c2xmin,c2xmax){ x=t; y=c2ymax ;label=wall;};
border center2l(t=c2ymin,c2ymax){ x=c2xmin; y=t ;label=wall;};

/*plot(bottom1(2*size)+bottom2(2*size)+bottom3(2*size)+bottom4(2*size)+bottom5(2*size)
    +rightb(2*size)+rightc(2*size)+rightt(2*size)
    +top1(2*size)+top2(2*size)+top3(2*size)+top4(2*size)+top5(2*size)
    +leftt(2*size)+lefth(2*size)+leftb(2*size)
    +centerb(size)+centerr(2*size)+centert(size)+centerl(2*size)
    +center2b(2*size)+center2r(size)+center2t(2*size)+center2l(size),wait=true);*/

//-------------- définition fonction flux ---------------------

func real fluxt(real t) { return (1.-(sin(abs(t-0.5*(rtmin+rtmax))/(rtmax-rtmin)*pi/2))^2);}
func real fluxb(real t) { return (1.-(sin(abs(t-0.5*(rbmin+rbmax))/(rbmax-rbmin)*pi/2))^2);}

//------------------ définition maillage ----------------
mesh Th = buildmesh(bottom1(2*size)+bottom2(2*size)+bottom3(2*size)+bottom4(2*size)+bottom5(2*size)
    +rightb(2*size)+rightc(2*size)+rightt(2*size)
    +top1(2*size)+top2(2*size)+top3(2*size)+top4(2*size)+top5(2*size)
    +leftt(2*size)+lefth(2*size)+leftb(2*size)
    +centerb(size)+centerr(2*size)+centert(size)+centerl(2*size)
    +center2b(2*size)+center2r(size)+center2t(2*size)+center2l(size));

/*plot(Th,wait=true);*/

//------------------ définition espace EF -------------
fespace Vh(Th,P1);

//------------- formulation variationnelle : potentiel ---------------

Vh phi, w;

solve potentiel(phi,w) =  int2d(Th)( dx(phi)*dx(w) + dy(phi)*dy(w) )
- int1d(Th,flowout)( paramOpt[2]*w )
+ on(flowin,phi=0);


// plot(phi,fill=true,value=true,wait=1);

//-------------- calcul des vitesse -----------------

Vh velx, vely;

velx = dx(phi);
vely = dy(phi);

//plot([velx,vely],wait=1);


//------------- formulation variationnelle : chaleur ---------------

Vh temp;

solve chaleur(temp,w) =  int2d(Th)( coeffdiff*(dx(temp)*dx(w) + dy(temp)*dy(w)) )
+ int2d(Th)( (velx*dx(temp) + vely*dy(temp) )*w)
- int1d(Th,radb)( coeffdiff*paramOpt[0]*fluxb(x)*w ) - int1d(Th,radt)( coeffdiff*paramOpt[1]*fluxt(x)*w )
 + on(hot,temp=temphot)
 + on(cold,temp=tempcold);

//plot(temp,fill=true,value=true,wait=true);

//------------- fonction coût -----------------

real fcost = 0.5*int2d(Th)( (temp - target )^2 );

cout  << "cost function= " << fcost << endl;


Vh phiAlpha;

solve potentielSensiPhiO(phiAlpha,w) =  int2d(Th)( dx(phiAlpha)*dx(w) + dy(phiAlpha)*dy(w) )
- int1d(Th,flowout)(w) /*paramOpt[2] c'est le potentiel*/
+ on(flowin,phiAlpha=0);


//plot(phiAlpha,fill=true,value=true,wait=1);

Vh velxAlpha, velyAlpha;

velxAlpha = dx(phiAlpha);
velyAlpha = dy(phiAlpha);

Vh tempAlpha1;

solve chaleurSensiPhiO(tempAlpha1,w) =  int2d(Th)( coeffdiff*(dx(tempAlpha1)*dx(w) + dy(tempAlpha1)*dy(w)) )
+ int2d(Th)( (velx*dx(tempAlpha1) + vely*dy(tempAlpha1))*w) +int2d(Th)((velxAlpha*dx(temp) + velyAlpha*dy(temp))*w  )
 + on(hot,tempAlpha1=0)
 + on(cold,tempAlpha1=0);

// plot(tempAlpha1,fill=true, value=true, wait=1); /* il montre que si on souffle plus au niveau du radiateur, la température va beaucoup diminuer, ça va entrainer l'air froid*/
/* la sensi sur phio va montrer comment la temperature evolue si on augmente le courant d'air*/

Vh tempAlpha2;

solve chaleurSensiPhiR1(tempAlpha2,w) =  int2d(Th)( coeffdiff*(dx(tempAlpha2)*dx(w) + dy(tempAlpha2)*dy(w)) )
+ int2d(Th)( (velx*dx(tempAlpha2) + vely*dy(tempAlpha2))*w) - int1d(Th, radt)( coeffdiff*fluxt(x)*w )
 + on(hot,tempAlpha2=0)
 + on(cold,tempAlpha2=0);

// plot(tempAlpha2,fill=true, value=true, wait=1);
 /*Si on augmente la temperature au niveau du radiateur 1 (en bas)
le graphe nous montre que la temperature va augmenter fort au niveau du radiateur et moins ailleurs, et ca reste proche de zero à cote des fenetres car 
valeurs fixées*/

Vh tempAlpha3;

solve chaleurSensiPhiR2(tempAlpha3,w) =  int2d(Th)( coeffdiff*(dx(tempAlpha3)*dx(w) + dy(tempAlpha3)*dy(w)) )
+ int2d(Th)( (velx*dx(tempAlpha3) + vely*dy(tempAlpha3))*w) - int1d(Th, radb)( coeffdiff*fluxb(x)*w )
 + on(hot,tempAlpha3=0)
 + on(cold,tempAlpha3=0);

//plot(tempAlpha3,fill=true, value=true, wait=1);
 /*c'est pas symetrique car il y a le courant d'air
en haut à droite*/

/* Calcul du gradient sensibilite*/

real gradJAlpha1 = int2d(Th)((temp - target )*tempAlpha1);
real gradJAlpha2 = int2d(Th)((temp - target )*tempAlpha2);
real gradJAlpha3 = int2d(Th)((temp - target )*tempAlpha3);

cout.precision(15);
cout << "Valeur fonction cout : " << fcost << endl;
cout << "Gradient de J par rapport a Alpha (sensibilite phi0) : " << gradJAlpha1 << endl; /* Sensibilite beaucoup + forte*/
cout << "Gradient de J par rapport a Alpha (sensibilite R1) : " << gradJAlpha2 << endl;
cout << "Gradient de J par rapport a Alpha (sensibilite R2) : " << gradJAlpha3 << endl;

/* Equations adjointes */

real lr = 1;  // learning rate à faire varier
int Nmax = 80; // nombre d'epochs à faire varier 
Vh tempopti;
Vh lambda1;
Vh lambda2;
real[int] funcCost(Nmax);
real alpha = 1.2;
real beta = 0.5;
real[int] it(Nmax);

for(int i=0; i<Nmax; i++){

solve potentiel(phi,w) =  int2d(Th)( dx(phi)*dx(w) + dy(phi)*dy(w) )
- int1d(Th,flowout)( paramOpt[2]*w )
+ on(flowin,phi=0);


// plot(phi,fill=true,value=true,wait=1);

//-------------- calcul des vitesse -----------------

velx = dx(phi);
vely = dy(phi);

// plot([velx,vely],wait=1);


//------------- formulation variationnelle : chaleur ---------------

solve chaleur(tempopti,w) =  int2d(Th)( coeffdiff*(dx(tempopti)*dx(w) + dy(tempopti)*dy(w)) )
+ int2d(Th)( (velx*dx(tempopti) + vely*dy(tempopti) )*w)
- int1d(Th,radb)( coeffdiff*paramOpt[0]*fluxb(x)*w ) - int1d(Th,radt)( coeffdiff*paramOpt[1]*fluxt(x)*w )
 + on(hot,tempopti=temphot)
 + on(cold,tempopti=tempcold);


solve eqAdjointe1(lambda2, w) = int2d(Th)((dx(phi)*dx(w) + dy(phi)*dy(w))*lambda2)
+ int2d(Th)(coeffdiff*(dx(w)*dx(lambda2) + dy(w)*dy(lambda2)))
+ int2d(Th)((tempopti - target ) * w)
+ on(hot,lambda2 = 0) /*condition aux limites (ici sur gammaF) comme precedemment*/
+ on(cold,lambda2 = 0);


solve eqAdjointe2(lambda1, w) = int2d(Th)((dx(tempopti)*dx(w) + dy(tempopti)*dy(w)) * lambda2) + int2d(Th)(dx(w)*dx(lambda1) + dy(w)*dy(lambda1))
+ on(flowin,lambda1 = 0); /*ici sur gammaI*/

/* Calcul du gradient (eq adjointe) */

real gradJAdjoint1 = -int1d(Th, flowout)(lambda1);
real gradJAdjoint2 = -int1d(Th, radb)(coeffdiff * fluxb(x) * lambda2);
real gradJAdjoint3 = -int1d(Th, radt)(coeffdiff * fluxt(x) * lambda2);
  
cout << "Gradient de J par rapport a Alpha (adjointe phi0) : " << gradJAdjoint1 << endl; /* Sensibilite beaucoup + forte*/
cout << "Gradient de J par rapport a Alpha (adjointe R1) : " << gradJAdjoint3 << endl;
cout << "Gradient de J par rapport a Alpha (adjointe R2) : " << gradJAdjoint2 << endl;

/* On a maintenant les valeurs de notre gradient, on peut donc faire de l'optimisation */
/* On privilégie l'eq adjointe car plus rapide et plus efficace */

/* Algorithme de descente du gradient */

real norm = sqrt(gradJAlpha2^2 + gradJAlpha1^2 + gradJAlpha3^2);
	
paramOpt[0] = paramOpt[0] - lr * gradJAlpha2/norm;
paramOpt[1] = paramOpt[1] - lr * gradJAlpha3/norm;
paramOpt[2] = paramOpt[2] - lr * gradJAlpha1/norm;

funcCost[i] = 0.5*int2d(Th)( (tempopti - target )^2 );

if(i != 0){
	if(funcCost[i] < funcCost[i-1]){
		lr = alpha * lr;
	}
	else{
		lr = beta * lr;
	}
}
  
}

for(int i=0; i<Nmax; i++){
	cout << i << " " << funcCost[i] << endl;
}

cout << "PhiR1 " << paramOpt[1] << endl;

cout << "PhiR2 " << paramOpt[0] << endl;

cout << "Phi0 : " << paramOpt[2] << endl;


plot(tempopti,fill=true, value=true, wait=1);
