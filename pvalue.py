import math as Math
def rprb(r,s,df):
    zp= [0,0,0]
    coef = [-0,1,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,1]
    EMAX=89
    idf=df
    dg=(df/2)-1
    idg=dg
    g=0
    xi=0
    dmax=0

    if (r==0.0):
        r=0.000000001
    if(dg <0):
        g=.572364943    
    elif (dg>0):
        if(dg== .5):
            thg=g-.120782238
        else:
            for i in range(1,int(idg)+1):
                g+=Math.log(dg-i+1)                
            
            
            if(df%2 !=0 ):
                g=g-.120782238 #4
        
    c = Math.exp(.693147181+Math.log(.196349541 *df)*df/2-g) # 5
    dt = Math.log(r+1)*100
    if(df>dt): dmax=dt
    else: dmax=df
    
    ri=5.2+r*.5-Math.log(dmax)-1/Math.sqrt(r)
    if (ri<.15):
        ri=.15
    bi= r-ri
    if(bi<0):
        bi=0
    ui=r+ri
    si=(ui-bi)/32
    sra=0
    for i in range(2,34):
        xi=i-1
        br = bi+si*xi
        sr =(-.918938533-br*br/(r*r *2)+Math.log(4*br/r))*df
        ra=0
        for j in range(1,34):
            xj=j-1
            bz=xj*.34375-7
            z = br+bz
            
            for k in range(1,3):
                if(k==2): 
                    z=bz
                x=abs(z)
                q = .39894228*Math.exp(-x*x/2)
                if(x>3.7): 
                    zp[k] = q*(Math.sqrt(4 +x*x)-x)/2
                
                
                if(x<=3.7):
                    t = 1./(1.+.2316419*x)
                    p=.31938153*t

                    p=p-.356563782*Math.pow(t,2)
                    
                    p=p+1.78147937*Math.pow(t,3)
                    p=p-1.821255978*Math.pow(t,4)
                    p=p+1.330274429*Math.pow(t,5)
                    
                    zp[k]=q*p
                
                xyz=0
                if(z>0): zp[k]=1-zp[k]
            
            rh=0
            d = zp[1]-zp[2]
            if(d>0):
                
                h = Math.log(d)
                pr = -.918938533-bz*bz/2 +h*(s-1)+Math.log(s)
                if(pr+EMAX>0): rh = Math.exp(pr)
                
        
            ra+=coef[j]*rh #10
        
    
        ra=.34375*ra/3
        srh=0
        d=(1-ra)/br
        if (d>0):
            sr+=Math.log(d)
            if((sr+EMAX)>0): srh=Math.exp(sr)
            
        sra = sra+coef[i]*srh #11
        
    result= si*sra/3*c
    return result