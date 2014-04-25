console.log('loaded USM align :-)');


align = function(s1,s2){
    if(!this.constructor.aligns){this.constructor.aligns={}}
    this.seq1 = s1;
    this.seq2 = s2;
    
    this.uid = 'align'+Math.random().toString().slice(2);
    this.constructor.aligns[this.uid]=this; // master list of align instances
    this.txt = function(){
        return ('\n'+s1+'\n'+s2+'\n');
    };
    
    this.showSeqs = function(){
        console.dir(this.seq1);
        console.dir(this.seq2);
        return [this.seq1.length,this.seq2.length]
    }
    
    // Identify alphabet    
    this.getABC = function(){
        var s12 = this.seq1+this.seq2;
        this.abc = jmat.unique(s12).sort();
    }
    
    // Smith-Waterman (http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
    // The simplified implementation, without the weighting of a score matrix
    this.smithWaterman=function(){
        console.log('Simplified Smith-Waterman');
        var t0=Date.now();
        this.getABC(); // get alphabet
        var n1 = this.seq1.length;
        var n2 = this.seq2.length;
        var s1 = ">"+this.seq1; // the initial ">" is just to force s[1] to be the first element ...
        var s2 = ">"+this.seq2; // ... such that the index 0 (zero) indicates the edge of H.
        var H = jmat.zeros(n1+1,n2+1); // scoring matrix
        var T = jmat.zeros(n1+1,n2+1); // traceback matrix
        var Hij = [0,0,0]; // initializaing matching score at position (i,j)
        var m = [0,0]; // initializing maximum, formated as [ max , index ], see jmat.max2 
        // Score
        for(var i=1;i<n1+1;i++){
            for(var j=1;j<n2+1;j++){
                Hij = [
                    H[i-1][j-1]+(s1[i]==s2[j]), // match
                    H[i-1][j], // insert in seq1
                    H[i][j-1]  // insert in seq2
                ];
                m = jmat.max2(Hij);
                H[i][j]=m[0];  // max -> score at position (i,j)
                T[i][j]=m[1];  // index -> to trace it back to the sequence
            }
        }
        // Trace it back
        var TB=[]; // traceback route
        var TBpos=[]; // position of traceback route
        i-=1;j-=1;
        while((i>0)&&(j>0)){
            TB.push(T[i][j]);
            TBpos.push([i,j]);
            //console.log(i,j," > ",T[i][j],H[i][j]);
            switch(T[i][j]){
                case 0:
                    i-=1;
                    j-=1;
                    break;
                case 1:
                    i-=1;
                    break;
                case 2:
                    j-=1;
                    break;
                default:
                    error('traceback not recognized !');
            }
        }
        TB.push(T[i][j]);
        TBpos.push([i,j]);
        //console.log(i,j," > ",T[i][j],H[i][j]);
        // finish traceback by sliding along the edge of H
        while(!((i==0)&&(j==0))){
            if(i>0){i-=1;TB.push(1)}
            else{j-=1;TB.push(2)}
            TBpos.push([i,j]);
        }
        TBpos.reverse(); TB.reverse(); // to get traceback results ordered forward
        //console.log("\n---------\n",i,j,"\n",JSON.stringify(H),"\n",JSON.stringify(T),"\n",JSON.stringify(TB));
        // Record results
        var score = H[H.length-1][H[0].length-1]; // total alignment score is last element of score matrix
        var t = Date.now()-t0;
        this.smithWaterman.H=H;
        this.smithWaterman.T=T;
        this.smithWaterman.score=score;
        this.smithWaterman.t=t;
        this.smithWaterman.TB=TB;
        this.smithWaterman.TBpos=TBpos;
    }
    
    
    
    // ---- UI stuff now, no science from here on ----
    
    
    this.UI=function(){ // create User Interface
        var div0 = document.getElementById('USMalign');
        if(!div0){
            div0 = document.createElement('div');
            div0.id = 'USMalign';
            div0.classList.add('container'); // in case bootstrap is available
            document.body.appendChild(div0);
        }
        var div = document.createElement('div');
        div.align=9; // pointer to the align object
        div.id = this.uid;
        div0.appendChild(div);
        div.appendChild(document.createElement('hr'));
        this.div=div; // it can always be found here
        div.textContent = 'Enter or paste 2 sequences: ';
        var demoSeqs = document.createElement('button');demoSeqs.id='demoSeqs';
        div.appendChild(demoSeqs);demoSeqs.textContent="demo seqs";demoSeqs.style.border=0;
        var in1 = document.createElement('textarea');in1.id="in1_"+div.id;
        in1.style.width='100%';in1.style.color='blue';
        var p1 = document.createElement('p');p1.appendChild(in1);div.appendChild(p1);
        var in2 = document.createElement('textarea');in2.id="in2_"+div.id;
        in2.style.width='100%';in2.style.color='blue';
        var p2 = document.createElement('p');p2.appendChild(in2);div.appendChild(p2);
        var smithAlign = document.createElement('button');smithAlign.id='smithAlign';
        div.appendChild(smithAlign);smithAlign.textContent="Simplified Smith-Waterman";
        var freeAlign = document.createElement('button');freeAlign.id='freeAlign';
        div.appendChild(freeAlign);freeAlign.textContent="Alignment-free";
        div.appendChild(document.createElement('hr'));
        // demo sequences
        var demo1 = 'AAABCBCBCBCBCDDDABBZZA';
        var demo2 = 'DDDBCBCBCBCBCAAAQQBZZ';
        // Actions
        demoSeqs.onclick=function(){
            //in1.value='TTGAAAGAAAAACAATTTTGGAATCGTATATTAGAATTTGCACAAGAAAGACTGACTCGATCCATGTATGA';
            //in2.value='AGAGTCGTCCGATTTTAACAGAACAATTTTGGAATCGTATATTAGAATTTGCACAAGAAAAGAAAA';
            in1.value=demo1;
            in2.value=demo2;
        }
        smithAlign.onclick=function(){
            var uid = this.parentElement.id;
            var aln = align.aligns[uid]; // this align
            var div = this.parentElement;
            if(in1.value.length==0){in1.value=demo1};
            if(in2.value.length==0){in2.value=demo2};
            aln.seq1 = in1.value;
            aln.seq2 = in2.value;
            if(!document.getElementById('SWcalc_'+uid)){
                var SWcalcDiv=document.createElement('div');
                SWcalcDiv.id = 'SWcalc_'+uid;
                div.appendChild(SWcalcDiv);
            } else {var SWcalcDiv = document.getElementById('SWcalc_'+uid)};
            SWcalcDiv.innerHTML = '<span style="color:red"> processing ...</span>';
            aln.smithWaterman();
            // show results
            var pre = document.createElement('pre');
            pre.textContent += '----- simplified Smith-Waterman alignment -----';
            var L1=""; // first line, seq 1
            var L2=""; // third line, seq 2
            var s1 = "-"+aln.seq1, s2 = "-"+aln.seq2; // blank insert character
            for (var i = 1 ; i< aln.smithWaterman.TBpos.length;i++){
                switch(aln.smithWaterman.TB[i]){
                    case 0:
                        L1+=s1[aln.smithWaterman.TBpos[i][0]];
                        L2+=s2[aln.smithWaterman.TBpos[i][1]];
                    break;
                    case 1:
                        L1+=s1[aln.smithWaterman.TBpos[i][0]];
                        L2+=s2[0];
                    break;
                    case 2:
                        L1+=s1[0];
                        L2+=s2[aln.smithWaterman.TBpos[i][1]];
                    break;
                }
            }
            pre.textContent +='\nscore: '+aln.smithWaterman.score+' time: '+aln.smithWaterman.t+'ms';
            pre.textContent +='\n'+L1+'\n'+L2;
            
            SWcalcDiv.innerHTML=""; // reset display
            SWcalcDiv.appendChild(pre);
            pre.style.color='green';
            
            
            4;
        }
        
    }
    
    
    
    // --- INI ---
    //this.showSeqs();
}

