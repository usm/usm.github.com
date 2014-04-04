console.log('loaded USM align :-)');

align = function(s1,s2){
    
    this.uid = 'UID'+Math.random().toString().slice(2);
    
    this.txt = function(){
        return ('\n'+s1+'\n'+s2+'\n');
    };
    
    this.seqs = function(){
        console.dir(s1);
        console.dir(s2);
    }
    
    
    // ---- UI stuff, no science from hereon ----
    
    
    this.UI=function(){ // create User Interface
        var div0 = document.getElementById('USMalign');
        if(!div0){
            div0 = document.createElement('div');
            div0.id = 'USMalign';
            div0.classList.add('container'); // in case bootstrap is available
            document.body.appendChild(div0);
        }
        var div = document.createElement('div');
        div.id = this.uid;
        div0.appendChild(div);
        div.appendChild(document.createElement('hr'));
        this.div=div; // it can always be found here
        var in1 = document.createElement('textarea');in1.id="in1";in1.style.width='100%'
        var p1 = document.createElement('p');p1.appendChild(in1);div.appendChild(p1);
        var in2 = document.createElement('textarea');in2.id="in2";in2.style.width='100%'
        var p2 = document.createElement('p');p2.appendChild(in2);div.appendChild(p2);
        div.appendChild(document.createElement('hr'));
        
    }
    
    this.seq1 = s1;
    this.seq2 = s2;
    
    
    
    
    // --- INI ---
    this.seqs();
}