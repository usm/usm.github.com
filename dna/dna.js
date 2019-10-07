console.log('dna.js loaded');

(function(){ // going anonymous for the sake of packaging
    const dna = function(){
        console.log('dna.js initiatized at '+Date())
        this.date=Date()
        //ini
    }
    //
    if(typeof(window)=="object"){
        window.dna=dna // drop it in the global scope 
    }
})()