console.log('entropy.js loaded')
// look for initial triggers at usm.entropy.ini

usm.entropy={
    plotArray(cgr){ // plot cgr array
        
    },
    dobra:function(u){
        return u
    },
    ui:function(){ // user interface
        var div = document.getElementById('usmEntropy')
        div.innerHTML='<div><span id="msg">type or paste sequence, or click </span><button id="demo">demo</button></div><textarea id="sequence" style="width:100%;height:30%"></textarea><button id="calc">Calculate entropy</button> '
        var calc = document.getElementById('calc')
        calc.onclick=function(){
            var seq = sequence.value // sequence in text area
            var u = new usm(seq)

            u = usm.entropy.dobra(u)
        }
        var demo = document.getElementById('demo')
        demo.onclick=function(){
            $.get('Es.seq.txt',function(x){
                x = x.split(/[\n\r]+/)
                console.log('read fastA file: '+x[0])
                msg.textContent='fastA file with '+x[0]+' '
                sequence.value=x[1]
                calc.click()
                
            })
            
            
        }
    },
    ini:function(){
        console.log('usm.entropy.ini started')
        this.ui() // build user interface
        4

    }
}

//usm.entropy.ini()

