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

usm = function (seq,abc,pack){ // Universal Sequence Map
	
	this.loadFasta=function(url,callback,abc){// load long sequences from a FastA file
		// for example
		// url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna' <-- big bacteria
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Acinetobacter_ADP1_uid61597/NC_005966.fna';
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Streptococcus_pneumoniae_R6_uid57859/NC_003098.fna'
		// url='ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Streptococcus_phage_2972_uid15254/NC_007019.fna' <-- phage
		// u = new usm;u.loadFasta(url,function(x){console.log(x.length)})
		// if default proxy doesn't work try this one: jmat.webrwUrl='http://webrw.no.de'
		jmat.disp('---- USMapping started ----');
		thisUsm = this;
		if (!!abc){this.abc=abc}
		//jmat.webrwUrl='http://webrw.no.de';
		jmat.get(url,function(x){ // encode sequences
			if (!!callback){callback(x)} // in case this function was called with a callback do that first
			// let's encode it now
			jmat.disp(x[0]); // fastA head identification, everything else should be a long sequence
			//x=x.splice(0,100); // <-- while debugging
			var n = x.length;
			thisUsm.seq = x.splice(1,n).reduce(function(a,b){return a+b}); // concatenate whole sequence
			x=''; // to free memory
			var nn = thisUsm.seq.length;
			jmat.disp('... found '+nn+' units divited in '+n+' segments,');
			thisUsm.encodeLong(); // the seq will be picked from this.seq
			//thisUsm.encode();
			}
		)
		jmat.disp('using proxy '+jmat.webrwUrl+' to get fastA file '+url+' ...');
		return 'using proxy '+jmat.webrwUrl+' to get fastA file '+url+' ...'
	}
	
	this.encodeLong = function(seq,abc,pack){ // encoding long sequences by writting directly to the usm instance
        if (!this.seq){throw ('Sequence not provided')}
        if (!this.abc){
			jmat.disp('find alphabet ...');
			this.abc=this.alphabet();	
		};
		jmat.disp('... alphabet: '+this.abc);
        this.cube=[];
		jmat.disp('packing USM space ...')
        this.str2cube(pack,true);
		jmat.disp('...',this.cube);
		jmat.disp('filling USM space ...');
        var m = this.cube.length, n = this.seq.length, i = 0;
        this.cgrForward = [];this.cgrBackward = [];
		//for(i=0;i<n;i++){this.cgrForward[i]=[[0]];this.cgrBackward[i]=[[0]]}
		for(i=0;i<m;i++){this.cgrForward[i]=[];this.cgrBackward[i]=[]}
        for (i=0;i<m;i++){
			jmat.disp((i*2+1)+'/'+(this.cube.length*2)+' mapping axis <'+this.cube[i]+'> forward');
            this.cgrLong(i,'cgrForward');
			jmat.disp((i*2+2)+'/'+(this.cube.length*2)+' mapping axis <'+this.cube[i]+'> backward');
			this.cgrLong(i,'cgrBackward');
            this.bin[i]=[]; // free memory
        }
		jmat.disp('packing CGR coordinates ...')
		jmat.disp('... forward ...')
		var c=[];
		for(i=0;i<n;i++){c[i]=[];for(j=0;j<m;j++){c[i][j]=this.cgrForward[j][i]}}
		this.cgrForward=c;
		jmat.disp('... backward ...')
		var c=[];
		for(i=0;i<n;i++){c[i]=[];for(j=0;j<m;j++){c[i][j]=this.cgrBackward[j][i]}}
		this.cgrBackward=c.reverse();
		var c=[];
		jmat.disp('---- USMapping done ----');
    }

	this.cgrLong=function(ii,direction){ // one dimension at a time
		var bin = this.bin[ii];
		var n = bin.length;//, bins=[], k=10; // binning threads
		if(direction=='cgrBackward'){bin.reverse();var c = this.cgrBackward[ii]}
		else{var c = this.cgrForward[ii]}
		var s=Math.random();for(var i=n-128;i<n;i++){s=s+(bin[i]-s)/2}
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

    this.transpose = function(M){
        var T=[];
        for (var i=0;i<M[0].length;i++){
            T[i]=[];
            for (var j=0;j<M.length;j++){T[i][j]=M[j][i]}
        }
        return T;
    }

	this.transpose2 = function(M){
        M[0].map(function(mi,i){var MM=[];})
    }

    this.encode = function(seq,abc,pack){
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
        for (var i=0;i<this.bin.length;i++){
            this.cgrForward[i]=this.cgr(this.bin[i]);
            this.cgrBackward[i]=this.cgr(this.bin[i].slice().reverse()).slice().reverse();
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

    if (seq){this.encode(seq,abc,pack)}
	else {if (abc){this.encode(abc,abc,pack)}} // if alphabet is provided then identify cube anyway to enable decoding

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
