// Universal Sequence Maps (hoding Object for CGR operations 

console.log('CGR toolbox :-)');

// Overloading String and Array - watch out!

String.prototype.reverse=function(){ // sort characters of a string
    return this.split('').reverse().toString().replace(/,/g,'');
}

String.prototype.sort=function(){ // sort characters of a string
    return this.split('').sort().toString().replace(/,/g,'');
}

Array.prototype.max=function(){ // returns maximum value
    return this.reduce(function(x1,x2){if(x1>x2){return x1}else{return x2}});
}

Array.prototype.min=function(){ // returns maximum value
    return this.reduce(function(x1,x2){if(x1<x2){return x1}else{return x2}});
}

Array.prototype.sum=function(){ // returns sum of all values
    return this.reduce(function(x1,x2){return x1+x2});
}

Array.prototype.transpose=function(){ // written for arrays of arrays
    if (!Array.isArray(this[0])){
        M=this.map(function(x){return [x]})}
    else {M=this};
    var T=[];
    for (var i=0;i<M[0].length;i++){
        T[i]=[];
        for (var j=0;j<M.length;j++){T[i][j]=M[j][i]}
    }
    if ((Array.isArray(T[0]))&&(T[0].length==1)){T=T.map(function(x){return x[0]})}
    return T;
}

// Universal Sequence Mapping (USM)

usm = function (seq,abc,pack,seed){ // Universal Sequence Map
	
	this.loadFasta=function(url,callback,abc){// load long sequences from a FastA file
		// for example
		// url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna' <-- big bacteria
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Acinetobacter_ADP1_uid61597/NC_005966.fna';
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Streptococcus_pneumoniae_R6_uid57859/NC_003098.fna'
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Streptococcus_phage_2972_uid15254/NC_007019.fna' <-- phage
		// u = new usm;u.loadFasta(url,function(x){console.log(x.length)})
		// if default proxy doesn't work try this one: jmat.webrwUrl='http://webrw.no.de'
		thisUsm = this;
		if (!!abc){this.abc=abc}
		//jmat.webrwUrl='http://webrw.no.de';
		//console.log('using proxy '+jmat.webrwUrl+' to get fastA file '+url);
		jmat.get(url,function(x){ // encode sequences
			if (!!callback){callback(x)} // in case this function was called with a callback do that first
			// let's encode it now
			console.log(x[0]); // fastA head identification, everything else should be a long sequence
			//x=x.splice(0,100); // <-- while debugging
			var n = x.length;
			thisUsm.seq = x.splice(1,n).reduce(function(a,b){return a+b}); // concatenate whole sequence
			x=''; // to free memory
			var nn = thisUsm.seq.length;
			console.log('... found '+nn+' units divited in '+n+' segments,');
			thisUsm.encodeLong(); // the seq will be picked from this.seq
			//thisUsm.encode();
			}
		)
		return 'using proxy '+jmat.webrwUrl+' to get fastA file '+url+' ...'
	}
	
	this.encodeLong = function(seq,abc,pack,seed){ // encoding long sequences by writting directly to the usm instance
		if(!!seq){this.seq=seq}
		if(!!abc){this.abc=abc}
        if (!this.seq){throw ('Sequence not provided')}
        if (!this.abc){
			console.log('find alphabet ...');
			this.abc=this.alphabet();	
		};
		console.log('... alphabet: '+this.abc);
        this.cube=[];
		console.log('packing USM space ...')
        this.str2cube(pack,true);
		console.log('...',this.cube);
		console.log('filling USM space ...');
        var m = this.cube.length, n = this.seq.length, i = 0;
        this.cgrForward = [];this.cgrBackward = [];
		//for(i=0;i<n;i++){this.cgrForward[i]=[[0]];this.cgrBackward[i]=[[0]]}
		for(i=0;i<m;i++){this.cgrForward[i]=[];this.cgrBackward[i]=[]}
        for (i=0;i<m;i++){
			console.log((i*2+1)+'/'+(this.cube.length*2)+' mapping axis <'+this.cube[i]+'> forward');
            this.cgrLong(i,'cgrForward');
			console.log((i*2+2)+'/'+(this.cube.length*2)+' mapping axis <'+this.cube[i]+'> backward');
			this.cgrLong(i,'cgrBackward');
            this.bin[i]=[]; // free memory
        }
		console.log('packing CGR coordinates ...')
		console.log('... forward ...')
		var c=[];
		for(i=0;i<n;i++){c[i]=[];for(j=0;j<m;j++){c[i][j]=this.cgrForward[j][i]}}
		this.cgrForward=c;
		console.log('... backward ...')
		var c=[];
		for(i=0;i<n;i++){c[i]=[];for(j=0;j<m;j++){c[i][j]=this.cgrBackward[j][i]}}
		this.cgrBackward=c.reverse();
		var c=[];
		console.log('USMapping done');
    }

	this.cgrLong=function(ii,direction,s){ // one dimension at a time
		var bin = this.bin[ii];
		var n = bin.length;//, bins=[], k=10; // binning threads
		if(direction=='cgrBackward'){bin.reverse();var c = this.cgrBackward[ii]}
		else{var c = this.cgrForward[ii]}
		if (!s){ // seed
			var s=Math.random();for(var i=n-128;i<n;i++){s=s+(bin[i]-s)/2}
		}
		//this[direction][ii]=[];
		//c[0][ii]=s+(bin[0]-s)/2;
		c[0]=s+(bin[0]-s)/2;
		for(var i=1;i<n;i++){
			//c[i][ii]=c[i-1][ii]+(bin[i]-c[i-1][ii])/2;
			c[i]=c[i-1]+(bin[i]-c[i-1])/2;
		}
		return 'done'
	}

    this.alphabet=function(seqString){//extracts alphabet
        if (!seqString){seqString=this.seq} //uses own string if argument not provided
        this.abc='';
        for (var i=0;i<seqString.length;i++){
            if (!this.abc.match(new RegExp(seqString[i]))){this.abc+=seqString[i]}
        }
        return this.abc.sort(); // using overloaded String.sort()
    }
    
    this.str2cube = function(pack,show){
        var m = this.abc.length;var n = this.seq.length;
        if (!pack){pack='compact'} // default packing method
        this.pack=pack;
        this.bin=[];
        switch (pack){
        case 'sparse':
            for (var j=0;j<m;j++){
                this.cube[j]=this.abc[j];
				this.bin[j]=this.seq.split('').map(function(si){
					if (si===this.abc[j]){return 0}
	                else {return 1}
				})
            }
            break;
        case 'compact':
            var L = Math.ceil(Math.log(m)/Math.log(2)); // map dimension
            var mm=Math.pow(2,L); // maximum length of this alphabet
            for (var j=0;j<L;j++){
                //this.bin[j]=[];
                var abc='';mm=mm/2;
                for (var i=0;i<m;i=i+mm*2){
                    abc+=this.abc.slice(i,i+mm);
                }
                this.cube[j]=abc;
				if (show){console.log((j+1)+'/'+L+' filling axis <'+abc+'>')}
				this.bin[j]=this.seq.split('').map(function(si){
				    if (abc.match(new RegExp(si))){return 0}
	                else {return 1}	
				})
            }
            break;
        //default :
        //  this.alphabet
        }
    }
        
    this.cgr = function(bin,s,ith){ // CGR with recursive seed
		if (!ith){ith=0}
        var n = bin.length;
		if (!s){s=bin[bin.length-1]} // start seed with last value of bin
        var y=[];
		y[0]=s+(bin[0]-s)/2;
		for(var i=1;i<n;i++){
			y[i]=y[i-1]+(bin[i]-y[i-1])/2;
		}
		// check recursive seed
		if ((s!=y[y.length-1])&(ith<64)){y=this.cgr(bin,y[y.length-1],ith+1)}
        return y;
    }

	this.cgr2 = function(bin,seed){ // CGR with 1/2 seed
		var y=[1/2+(bin[0]-1/2)/2];
		for(var i=1;i<bin.length;i++){y[i]=y[i-1]+(bin[i]-y[i-1])/2}
		return y;
	}

    this.transpose = function(M){
        var T=[];
		//try{
		for (var i=0;i<M[0].length;i++){
            T[i]=[];
            for (var j=0;j<M.length;j++){T[i][j]=M[j][i]}
        }
		//}
		//catch(err){
		//console.log(err);
		//  lala = 4;
		//}
        return T;
    }

	this.transpose2 = function(M){
        M[0].map(function(mi,i){var MM=[];})
    }

    this.encode = function(seq,abc,pack,seed){
        if (!seq){seq=this.seq}
        else {this.seq=seq}
        if (!this.seq){throw ('Sequence not provided')}
        if (abc){if(abc.length==0){var abc = undefined}} // in case abc=''
        if (abc){this.abc=abc}
        //if (!this.abc){this.abc=''}
        if (!this.abc){this.abc=this.alphabet()};
        var m = this.abc.length;var n = this.seq.length;
        this.cube=[];
        this.str2cube(pack);
        this.cgrForward = [];
        this.cgrBackward = [];
		if (!seed){
        	for (var i=0;i<this.bin.length;i++){
            	this.cgrForward[i]=this.cgr(this.bin[i]);
            	this.cgrBackward[i]=this.cgr(this.bin[i].slice().reverse()).slice().reverse();
        	}
		}
		else{ // seeded CGR
			for (var i=0;i<this.bin.length;i++){
            	this.cgrForward[i]=this.cgr2(this.bin[i],seed);
            	this.cgrBackward[i]=this.cgr2(this.bin[i].slice().reverse(),seed).slice().reverse();
        	}
		}
        //delete this.bin; // comment out if .bin is of no use <---<---<--- NOTE!
        this.cgrForward=this.transpose(this.cgrForward);
        this.cgrBackward=this.transpose(this.cgrBackward);
		this.usm=[];
		for (var i=0;i<n;i++){ // In a serious application .usm is all we'd need to keep .cgr--- etc could all be deleted
			this.usm[i]=[this.cgrForward[i],this.cgrBackward[i]];
		}
    }

    this.decodeBin = function(x,n){// decompose single coordinate, x, into a binary array of length <= n
        if (!n){n=Infinity}
        var y = [];var i = 0;
        x=x*2;var z=[];
        y[0]=x; // just in case x is 0, 1 or even 1/2 (impossible) y is still populated
        while ((x!=Math.round(x))&(i<n)){
            if (x>1){y[i]=1;x=x-1}
            else{y[i]=0}
            x = x*2;i++;
        }
        return y;
    }

    this.bin2int = function (B){//converts binary vector into integer
        return B.slice().reverse().map(function(x,i){return x*Math.pow(2,i)}).reduce(function(a,b){return a+b});
    }

    this.decode = function(xy,n){// decode numerical coordinates
        var bin2int = this.bin2int;
        var decodeBin = this.decodeBin;
        var abc = this.abc;
        var bb = xy.map(function(x){return decodeBin(x)});
        bb = this.transpose(bb).map(function(x){return bin2int(x)});
        bb = bb.map(function(x){return abc[x]});
        return bb.toString().replace(/,/g,'');
    }

    // run USM map automatically is a sequence is provided
    
    this.L = function (a,b){ // distance between two coordinates
        var d=0;
        while((Math.pow(2,d)!=Infinity)&Math.round(a*Math.pow(2,d))==Math.round(b*Math.pow(2,d))){d++}
		if(Math.pow(2,d)==Infinity){d=64} // stop at the numerical resolution
        return d;
    }

    this.distCGR = function (a,b){ // distance between two CGR positions, a and b
        //var ab=this.transpose([a,b]); // such that each position is an array of elements with two coordinates, one from each sequence
        var dist = this.L;
        return this.transpose([a,b]).map(function(x){return dist(x[0],x[1])}).min();
    }
    
    this.mapReduce = function (x,map,reduce){
        return reduce(x.map(map))
    }

    this.distProfile= function (s){// Distance to another sequence, s
        if (typeof(s)=='string'){
            s=new usm(s,this.abc,this.pack);
        }
        if (s.abc!==this.abc){throw('unequal alphabets')}
        if (s.pack!==this.pack){throw('unequal encode packing')}
        //var dist = this.distCGR;
        //var distCoord
        //var dist = function(a,b){return this.transpose([a,b]).map(function(x){return this.L(x[0],x[1])}).min()};
        var f = this.cgrForward.map(function(x){var x0=x; return s.cgrForward.map(function(x){return s.distCGR(x0,x)}).sum()});
        var b = this.cgrBackward.map(function(x){var x0=x; return s.cgrBackward.map(function(x){return s.distCGR(x0,x)}).sum()});

        return [f,b].transpose().sum();
    }

	this.dist = function (x,y){ // Equation 5: distance between two usm coordinates for a position, i
		                        // each provided as a two element Array [cgrForward[i],cgrBackward[i]]
		var d=this.distCGR(x[0],y[0])+this.distCGR(x[1],y[1])-1;
		if (d<0){d=0}
		return d
	}

	this.distMap=function (sprobe,sbase){// Distance Map to a new probing sequence, sprobe
		if (typeof(sprobe)=='string'){sprobe = new usm(sprobe,this.abc,this.pack)} // in case the actual sequence was inputed
		if (!sbase){sbase=this} //
		return sbase.usm.map(function(x){return sprobe.usm.map(function(y){return sbase.dist(x,y)})});
	}
	
	this.align=function(sprobe,sbase){ //align sequence sprobe to this sequence
		if (typeof(sprobe)=='string'){sprobe = new usm(sprobe,this.abc,this.pack)} // in case the actual sequence was inputed
		if (!sbase){sbase=this} //
		// start by considering a complete match and zoom in until such a subset is found
		var n = sprobe.seq.length , nn = sbase.seq.length;
		var A={posBase:[0],posProbe:[0],match:[0],ind:0}
		this.alignUsm=function(posStart,posEnd){ // defined here to capture align closure
			console.log([posStart,posEnd]);
			var inc = Math.floor((posEnd-posStart)/2);// increment
			var i=posStart+inc;
			var res=32;//some resolution
			//this.distCGR(x[0],y[0])+this.distCGR(x[1],y[1]) Eq 5
			var d = sbase.cgrForward.map(function(cf,ii){return sbase.distCGR(cf,sprobe.cgrForward[i])+sbase.distCGR(sbase.cgrBackward[ii],sprobe.cgrBackward[i])});
			var mm = jmat.max2(d); // lets code it here so we don't need jmat just yet:
			//var mm = d.reduce(function(a,b){return Math.max(a,b)});
			if(mm[0]>0){mm[0]-=1};// var ind = mm[0],dmax=mm[1];
			var j=A.posProbe.length;
			//var dF=sbase.distCGR(sbase.usm[mm[1]][0],sprobe.usm[i][0]); // forward distance
			var dF=sbase.distCGR(sbase.cgrForward[mm[1]],sprobe.cgrForward[i]); // forward distance
			if(dF>0){dF-=1}
			A.posBase[j]=mm[1]-dF; // using foward CGR to find start position
			if(A.posBase[j]<0){mm[0]-=(dF-mm[1]);dF=mm[1];A.posBase[j]=0} // check lower boundary for Base
			if((A.posBase[j]+dF)>nn){mm[0]-=A.posBase[j]+mm[0]-nn} // check upper boundary for Base
			A.posProbe[j]=i-dF;
			if(A.posProbe[j]<0){mm[0]-=(dF-i);dF=i;A.posProbe[j]=0;A.posBase[j]=mm[1]-dF;} // check lower boundary for Probe
			if((A.posProbe[j]+mm[0])>n){mm[0]-=A.posProbe[j]+mm[0]-n} // check upper boundary of Probe
			if(mm[0]>res){ // check if numerial resolution might have been exceeded
							var dFi=sbase.distCGR(sbase.cgrForward[A.posBase[j]],sprobe.cgrForward[A.posProbe[j]])-1; // additional forward distance
							while(dFi>0){ // if some was found
								A.posProbe[j]-=dFi;A.posBase[j]-=dFi;
								mm[0]+=dFi;
								dFi=sbase.distCGR(sbase.cgrForward[A.posBase[j]],sprobe.cgrForward[A.posProbe[j]])-1;
							}
							dFi=sbase.distCGR(sbase.cgrBackward[A.posBase[j]+mm[0]-1],sprobe.cgrBackward[A.posProbe[j]+mm[0]-1])-1; // backward distance at the end of match segment
							while(dFi>0){ // if some was found
								mm[0]+=dFi;
								dFi=sbase.distCGR(sbase.cgrBackward[A.posBase[j]+mm[0]-1],sprobe.cgrBackward[A.posProbe[j]+mm[0]-1])-1;
							}
						}
			A.match[j]=mm[0];
			if(mm[0]>(A.match[A.ind])){A.ind=j}; // if this is the best match yet
			if(A.match[A.ind]<inc){ // if best match is shorter than increment, zoom into its halves
				inc = Math.floor(inc/2);
				this.alignUsm(posStart,inc);
				this.alignUsm(inc+1,posEnd);
			}
		}
		this.alignUsm(0,n);
		console.log('largest identical segment has length '+A.match[A.ind]+' and aligns with position '+A.posBase[A.ind]+' in base sequence and position '+A.posProbe[A.ind]+' in probe sequence');
		return A
		//sbase.usm.map(function(x){return sprobe.usm.map(function(y){return sbase.dist(x,y)})});
	}

    this.alignQ = function (sprobe, sbase) {
     // This function provides an asynchronous wrapper for the usual
     // `this.align` method using Quanah (http://wilkinson.github.com/quanah/).
        if (Object.prototype.hasOwnProperty('Q') === false) {
            throw new Error('Quanah is not loaded.');
        }
        var Q, that, y;
        Q = Object.prototype.Q;
        that = this;                    //- the current USM object
        y = Q.avar();
        y.Q(function (evt) {
            if (that.hasOwnProperty('cgrBackward') === false) {
                window.setTimeout(y.revive, 0);
                return evt.stay('Waiting for indexing to finish ...');
            }
            y.val = that.align(sprobe, sbase);
            return evt.exit();
        });
        return y;
    };

    if (seq){// find out if this is a sequence or the url of a fastA file with one
		if(seq.length>15){ // it could be a url
			if(!!seq.slice(0,10).match(/:\/\//)){
				console.log('using proxy '+jmat.webrwUrl+'\nto get fastA file '+seq+' ...');
				this.loadFasta(seq,function(x){console.log('... file loaded ...')}, abc)}
			else{this.encode(seq,abc,pack)}
		}
		else{this.encode(seq,abc,pack,seed)}
	}
	else {if (abc){this.encode(abc,abc,pack,seed)}} // if alphabet is provided then identify cube anyway to enable decoding

}

usmDist=function(s1,s2,opt){ // Compares two USM encoded sequences
	if (!opt){opt='matrix'};
	switch (opt)
	{
	case 'matrixForward':
		return s1.cgrForward.map(function(x){var x0=x; return s2.cgrForward.map(function(x){return s2.distCGR(x0,x)})});
		break;

	case 'matrixBackward':
		return s1.cgrBackward.map(function(x){var x0=x; return s2.cgrBackward.map(function(x){return s2.distCGR(x0,x)})});
		break;

	case 'matrix':
		//var f1,b1,f2,b2;
		//create two mirror (forward and backward) version of s1 and s2
		fb1=s1.transpose([s1.cgrForward,s1.cgrBackward]);
		fb2=s2.transpose([s2.cgrForward,s2.cgrBackward]);
		return fb1.map(function(x){var x0=x; return fb2.map(function(x){var res=s1.distCGR(x0[0],x[0])+s1.distCGR(x0[1],x[1]);if(res>1){res=res-1};return res})});
		break;

	case 'sum':
		break;

	case 'max':
		break;

	case 'min':
		break;


	}
}

// load jmat
if(typeof(jmat)=='undefined'){
	var s=document.createElement('script');
	s.src='https://jmat.googlecode.com/git/jmat.js';
	document.head.appendChild(s);
}
