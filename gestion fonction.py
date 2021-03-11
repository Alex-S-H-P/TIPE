import random
import matplotlib.pyplot as plt
factor = 1.7


def creer(a = 0 ,b = 0,f_a = None,f_b = None):
    """renvoie une fonction sous forme dict 
    la fonction définie entre a et b qui renvoie f_a et f_b à ses bornes respecives."""
    if a > b:
        a  ,  b  = b,a
        f_a, f_b = f_b, f_a
    
    if f_a == None:
        f_a = random.normalvariate(0,10**factor)
    if f_b == None:
        f_b = random.normalvariate(0,10**factor)
    
    dic = {a:f_a,b:f_b}
    return dic

def fonction(dic):
    """ prends un dictionnaire de valeurs, et renvoie f"""
    f  = 'fonction_cache_124578963 = lambda x : '
    for a in dic:
        value = dic[a]
        f += str(value) + '*('
        g = ''
        h = ''
        for b in dic:
            if b != a:
                g += ' (x-' +     str(b) +      ') *'
                h += ' (' + str(a)+'-'+ str(b)+' ) *'
        g = g[:-1]
        h = h[:-1]
        f += g + ')/(' + h + ') +'
    f = f[:-1]
    print(f[27:])
    exec( 'global fonction_cache_124578963 \n' + f )
    return fonction_cache_124578963
    
def iteration(sample):
    ## establishing the cost
    cost_list = []
    for f in sample:
        cost_list.append(cost(f))
    ## rangement de la chambre
    sample , cost_list = ranger(sample,cost_list)
    n = len(sample)
    i = n // 2
    ## meurtre
    sample = sample [:i]
    ## reproduction
    for j in range(i):
        sample.append(modify(sample[i]))
    while len(sample) < n:
        sample.append(generate())
    return sample,list_cost

def algo_gen(n,score_lim):
    ## initialisation
    sample = []
    for i in range(n):
        sample.append(generate())
    Y1,Y2,X,x,score = [],[],[],0,[0]
    ## iteration
    while score[0] < score_lim:
        Y1.append( score[0] )
        Y2.append(moy(score))
        X.append(x)
        sample,score = iteration(sample)
        x += 1
        ## effets visuels
        plt.plot(X,Y1)
        plt.plot(X,Y2)
        plt.show()
    ## fin
    return sample[0]
    
def moy(iterable):
    """returns an avarged number out of the iterable"""
    sum = 0
    for truc in iterable:
        sum += truc
    return sum / len(iterable)

def ranger(liste, valeur):
    """renvoie le couple/liste/,/valeur/ rangés pour valeur décroissante"""
    for i in range(len(liste)):
        item = {"value" : valeur [i], " data " : liste[i]
        for j in range(i,len(liste)):
            if valeur[j] > item["value"]:
                item["value"],valeur[j] = valeur[j],item["value"]
                item["data" ],liste [j] = liste [j],item["data" ]
        valeur[i] = item["value"]
        liste [i] = item["data" ]
    return liste,valeur

def generate():
    """ renvoie un objet légal """
    
                
