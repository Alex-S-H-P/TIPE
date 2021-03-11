import matplotlib.pyplot as plt
import numpy as np
import time as T
from threading import Thread
pi = np.pi
g =  9.80665

a = .5
w = np.pi 
waveheight = .5
epsilon = 1/12          # epsilon correspondra, tout au long du programme, à la précision.
    
yes = ('y','ready','Yes','o','oui','Oui','yes')
no  = ('n','nope' ,'No','no','non','Non','Nope')

def entre(valeur,mini,maxi):
    """borne la /valeur/ entre /min/ et /max/"""
    return min(max(valeur,mini),maxi)

def dist(A,B):
    """renvoie la distance entre A = (xA,yA) et B = (xB, yB)"""
    return ((B[0]-A[0])**2 + (B[1]-A[1])**2 )**(.5)

def integr( f , borne_i , borne_f , e = epsilon):
    """réalise l'intégrale d'une fonction /f/ entre /borne_i/ et /borne_s/ avec un écart e de base 10**-8"""
    try:
        t = borne_i
        I = 0
        nb_points = 0
        while t < borne_f :
            t += e/2
            I += f(t)
            t += e/2
            nb_points += 1
        return I*e
    except:
        print(f)
        return 0

def f(O , X , t):
    x = dist(O,X)
    return waveheight * (np.cos(-a * x + w * -t) + 0.5 * np.cos(2 * (-a * x )+ w * -t))

def signe(valeur):
    """renvoie 1 si valeur > 0, 0 si valeur = 0 et -1 sinon"""
    
    if valeur > 0:
        return 1
    elif valeur == 0:
        return 0
    elif valeur < 0:
        return -1
    else:
        reason = 'unexpected value : value does not have defined sign : ' + str(valeur)
        raise invalide(reason)

def Id(element,liste):
    assert element in liste, "{} n'est pas dans la liste {}".format(element,liste)
    for i in range(len(liste)):
        truc = liste[i]
        if truc == element:
            return i
        

def reverse(liste):
    """renvoie la liste à l'evers:"""
    l = []
    for e in liste:
        l = [e] + l
    return l

class invalide(Exception):
    
    def __init__(self , reason = None):
        self.message = 'Something went wrong'
        self.reason = reason
    def __str__(self):
        if self.reason is not None:
            return '{0} : The identified reason is {1}'.format(self.message, self.reason)
        else:
            return self.message

class world:
    
    def __init__(self, ciel = (116/255, 208/255, 241/255), mer = True):
        """
        Exemple de syntaxe de création d'objet: >>c2 = world.cube(1,monde)
        """
        self.t = 0
        self.ciel = ciel
        self.objets = []
        if mer:
            self.sea = world.mer(f)
            self.objets.append(self.sea)
    
    def update(self, dt):
        """met tout le monde à jour des modifications réalisées en /dt/ secondes.
        
        De base, merci de mettre tous les objets en objet"""
        for objet in self.objets:
            if not objet.__class__.__name__ in ('mer', 'sea'):
                objet.update(dt)
        self.t += dt
    
    class mer:
        def __init__(self, f , O = [-1000,0,0]):
            """f dépend de O , X et de t( dans cet ordre)
            O est le centre des vagues.
            X est un point du plan z = 0
            t est en secondes"""
            self.f = f
            self.O = O
        
        def F(self,X,t):
            return self.f(self.O,X,t)
        
        def niveau(self, cote, t):
            """renvoie le niveau moyen de l'eau dans une colone de section horizontale carrée
            
            /!\ Calcul de 1/(epsilon²) points
            """
            Σ,i = 0,0
            for k in range(int(1/epsilon)):
                for j in range(int(1/epsilon)):
                    y,x = (k*epsilon - 1/2)* cote, (j*epsilon - 1/2)* cote
                    Σ += self.F((x,y),t)
                    i += 1
            return Σ / i
        
        
    class Objet:
        
        def __init__(self, monde):
            ##  référencement
            monde.objets.append(self)
            self.monde = monde
            ## inititialisation positonnelle
            self.z = 0
            self.v_z = 0
            self.P = 0
            self.Π = 0
            ## jusqu'à nouvel ordre
            self.Cx = 1.2
            self.result = None
        
        def __str__(self):
            """ce qui sera affiché lors de la converion de l'objet en str"""
            return 'Object : ' + self.__class__.__name__ + '| Stored as : ' + str(self.__class__.__mro__) + ' in ' + str(self.__sizeof__()) + ' bytes'
            
        def O(self):
            return 0,0
        
        def update(self , dt , frottements_fluides = False):
            """réalise l'opération //objet//(t) -> //objet//(t+dt)
            
            le boolean frottements_fluides fait appel à la méthode du même nom de l'objet (non définit ici).
            Attention à mettre la simulation à jour avant de lancer cette page de l'histoire.
            
            On réalise ici:
                            dv = Σ(F) * dt
                            dz =  dv  * dt
            
            """
            z_e = self.monde.sea.niveau(self.cote,self.monde.t)
            self.Π = self.archimede(z_eau = z_e)
            self.P = self.poids()
            ΣF = self.Π - self.P
            if frottements_fluides:
                ΣF += self.FF(z_e)
            self.v_z += (ΣF * dt)/self.masse
            self.z += self.v_z * dt
            if abs(self.z) >= 50:
                self.z = signe(self.z)*50
            if abs(self.v_z) >= 2**20:
                print(self.v_z)
                self.v_z = 2 ** 20
        
        def _play_(self,ti,tf):
            """Réalise la simulation de l'objet dans le monde"""
            ## simulation
            self.z = 0
            self.v = 0
            reponse = {"min":0,"max":0,"v_max":0}
            self.monde.t = 0
            T,Z = [],[]
            while self.monde.t <= tf:
                self.monde.update(epsilon)
                if self.monde.t >= ti:
                    Z.append(self.z)
                    T.append(self.monde.t)
                if self.z > reponse["max"]:
                    reponse["max"] = self.z
                elif self.z < reponse["min"]:
                    reponse["min"] = self.z
                if abs(self.v_z) > reponse["v_max"]:
                    reponse["v_max"] = abs(self.v_z)
            ## enregistrement
            self.result = reponse
        
        def amplitude(self,ti, tf):
            """Renvoie l'amplitude du mouvement de l'objet avec le monde dans les conditiuons que celui_ci a."""
            if self.result == None:
                self._play_(ti,tf)
            result = self.result
            return abs(result["max"]-result["min"])
            
        
        def v_max(self,ti,tf):
            """Renvoie la vitesse max (en valeur absolue) subie par l'objet"""
            if self.result == None:
                self._play_(ti,tf)
            result = self.result
            return result["v_max"]
                
        
        def FF(self,z_e,p_eau = 997):
            """réalise le calcul des Frottements Fluides (d'où FF)
            
            On utilise ici la formule quadratique de la vitess(plutôt appropriée pour les vitesses moyennes à importantes
            """
            niveau = z_e - (self.z-self.cote)/2
            niveau = entre ( niveau , 0, self.cote) # si le niveau est en dessous, ou au dessus, on prends juste la valeur limite
            V_em = 0
            return -signe(self.v_z) * 1/2 * self.Cx * (self.v_z )** 2 * self.section(niveau) *  p_eau
        
        def archimede(self , p_fluide = 997 , z_eau = 0 ):
            """calcule la poussée d'archimède du cube dans un fluide (de base, l'eau) de niveau /z_eau/, alors que l'objet a pour niveau /z_g/ (à son barycentre)"""
            if self.z >= z_eau:
                V_imm = 0
            else:
                V_imm = integr(self.section,self.z,min(z_eau,self.z + self.hauteur)) # le min prends en compte le cas où l'objet est immergé
            return V_imm * p_fluide * g
        
        def poids(self):
            global g
            return g*self.masse
        
        def set_Mv_to(self,masse_volumique):
            self.masse_volumique = masse_volumique
            self.masse = self.masse_volumique * self.volume
        
        
    class cube(Objet):
        
        def __init__(self, cote, monde, masse = 1):
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## dimensions
            self.cote = cote
            self.masse = masse
            self.volume = cote**3
            self.surface = 6*cote**2
            self.masse_volumique = self.masse / self.volume
            self.hauteur = self.cote
        
        def section(self, niveau):
            """renvoie la section du cube a un plan horizontal de hauteur dans le monde niveau"""
            niveau -= self.z
            if niveau < self.cote and niveau > 0 :
                return self.cote **2
            else:
                return 0
        
    class table(Objet):
        
        def __init__(self , dimensions , monde , masse = 1):
            """/self/ est une table de plateau carré et aux quatres coins des quels d'éventuels pieds (par défault non) sont disposés
            
            /dimensions/ est soit composé de:
                > largeur_plateau, hauteur_plateau . . . . . . . . . . . . . . . pour une table sans pieds
                > largeur_plateau, hauteur_plateau , cote_pied , hauteur_pied - -pour une table avec pieds
                
            /world/ correspond au monde d'appel, dont on pourra récupérer le temps /world.t/ ultérieurement.
            De plus, /world/ est munis d'une liste de tous les objets considérés.
            """
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## dimensions
            if len(dimensions) == 2:
                largeur_plateau, hauteur_plateau = dimensions
                cote_pied = 1
                hauteur_pied = 0 
            elif len(dimensions) == 4:
                 largeur_plateau, hauteur_plateau , cote_pied , hauteur_pied = dimensions
            else:
                error = invalide()
                raise invalide
                quit()
            ## attribution des valeurs
            cote_plateau = largeur_plateau
            self.cote_plateau = cote_plateau
            self.cote = hauteur_plateau + hauteur_pied
            self.cote_pied = cote_pied
            self.hauteur_pied = hauteur_pied
            self.largeur_plateau = largeur_plateau
            self.hauteur_plateau = hauteur_plateau
            self.surface = 16 * cote_pied * hauteur_pied + 4 * hauteur_plateau * cote_plateau + 2* (cote_plateau**2) #simplification par décalage de la surface des bas des pieds sur le bas du plateau (cf brouillon)
            self.volume = 4 * (cote_pied ** 2) * hauteur_pied + (cote_plateau ** 2) * hauteur_plateau
            self.masse = masse
            self.masse_volumique = self.masse / self.volume
            self.hauteur = self.hauteur_pied + self.hauteur_plateau
        
        def section(self,niveau):
            """renvoie la section du cube a un plan horizontal de hauteur dans le monde niveau"""
            niveau -= self.z
            if niveau < 0:
                return 0
            elif niveau <= self.hauteur_pied:
                return 4*self.cote_pied**2
            elif niveau  <= self.hauteur:
                return self.cote_plateau
            else:
                return 0
    
    class pillier(Objet):
        
        def __init__(self, dim, monde , masse = 1):
            """un pilier est un cube de coté /a/ pour lequel une bande de hauteur /a-2b/ de profondeur /c/ a été découpé sur 4 faces adjacentes (cf cahier)
             a , b , c = /dim/
            /world/ correspond au monde d'appel, dont on pourra récupérer le temps /world.t/ ultérieurement.
            De plus, /world/ est munis d'une liste de tous les objets considérés."""
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## dimensions
            (a,b,c) = dim
            self.masse = masse
            self.a = a
            self.cote = a
            self.hauteur = self.cote
            self.b = b
            self.c = c
            self.volume = (a**2) * (2*b)  + ((a-(2*c)) ** 2) * (a-(2*b))
            self.masse_volumique = self.masse / self.volume

        def section(self,niveau):
            """renvoie la section du cube a un plan horizontal de hauteur dans le monde niveau"""
            niveau -= self.z
            if (0 <= niveau and niveau <= self.b ) or (self.a-self.b <= niveau and niveau <= self.a):
                return self.a**2
            elif niveau >= self.b and self.a-self.b >= niveau:
                return (self.a-2*self.c)**2
            else:
                return 0
        
            
    class pillier2(Objet):
        
        def __init__(self, dim , monde , masse = 1):
            """
            /dim/ = a,b,c
            
            un pillier de catégorie 2 est un cube de coté /a/ pour lequel un creux carré de hauteur /a-2b/ de profondeur /c/ a été découpé sur 4 faces adjacentes (cf cahier)
            
            /world/ correspond au monde d'appel, dont on pourra récupérer le temps /world.t/ ultérieurement.
            De plus, /world/ est munis d'une liste de tous les objets considérés."""
            a,b,c = dim
            
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## dimensions
            self.masse = masse
            self.a = a
            self.cote = a
            self.hauteur = self.cote
            self.b = b
            self.c = c
            self.volume = (a**2) * (2*b)  + (((a-(2*c)) ** 2 ) + 4 * (c ** 2)) * (a-(2*b)) 
            self.masse_volumique = self.masse / self.volume

        def section(self,niveau):
            """renvoie la section du cube a un plan horizonbtal de hauteur dans le monde niveau"""
            niveau -= self.z
            if (0 <= niveau and niveau <= self.b ) or (self.a-self.b <= niveau and niveau <= self.a):
                return self.a**2
            elif niveau >= self.b and self.a-self.b >= niveau:
                return (self.a-2*self.c)**2 + 4 * (self.c**2)
            else:
                return 0

    class pillier3(Objet):
        
        def __init__(self, dim , monde, masse = 1):
            """
            (a,b,c, largeur_aile, longueur_aile)= /dim/
            
            crée un pillier
                        
            /world/ correspond au monde d'appel, dont on pourra récupérer le temps /world.t/ ultérieurement.
            De plus, /world/ est munis d'une liste de tous les objets considérés."""
            (a,b,c, largeur_aile, longueur_aile)= dim
            
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## dimensions
            self.masse = masse
            self.a = a
            self.cote = a
            self.hauteur = a
            self.b = b
            self.c = c
            self.l_aile = largeur_aile
            self.L_aile = longueur_aile
            self.section_aile = self.l_aile * self.L_aile
            self.volume = (a**2) * (2*b)  + ((a-(2*c)) ** 2 + 4 * self.section_aile) * (a-(2*b))
            self.masse_volumique = self.masse / self.volume

        def section(self,niveau):
            """renvoie la section du cube a un plan horizonbtal de hauteur dans le monde niveau"""
            niveau -= self.z
            if (0 <= niveau and niveau <= self.b ) or (self.a-self.b <= niveau and niveau <= self.a):
                return self.a**2
            elif niveau >= self.b and self.a-self.b >= niveau:
                return (self.a-2*self.c)**2 + 4 * self.section_aile
            else:
                return 0
        
    class defineable(Objet):
        
        def __init__(self,monde,f,niv_min = 0,niv_max = 20):
            ## initialisation systématique
            world.Objet.__init__(self,monde)
            ## initialistaion particulière
            g = fonction(f)
            self.volume  = integr(g,niv_min,niv_max)
            self.hauteur = abs(niv_max - niv_min)
            self.masse = self.volume
            ## fonction
            self.function = f
            self.S        = g
            ## paramètres
            self.min = niv_min
            self.max = niv_max
            self.set_masse()
            self.cote = self.S(niv_max)
        
        def clone(self):
            """renvoie un clone non lié de l'objet"""
            return world.defineable(self.monde,self.function,self.min,self.max)
        
        def mutate(self,f):
            """enregistre une nouvelle fonction /S(z)/"""
            self.function = f
            self.S        = fonction(self.function)
            self.set_masse()
            self.result = None
            return self
        
        def __str__(self):
            return str(self.function)
            
        def section(self,niveau):
            if   niveau < self.min:
                return 0 
            elif niveau > self.max:
                return 0
            else:
                try:
                    return max(self.S(niveau),0)
                except:
                    return 0
        
        def set_masse(self):
            function = self.section
            l , ρ = 0.2 , 1.2 *997
            g = lambda z : 4 * l * epsilon * (np.pi*function(z))**(.5)
            self.masse = integr(g,self.min,self.max,) * ρ
        
        def archimede(self , p_fluide = 997 , z_eau = 0 ):
            if self.z - self.min >= z_eau:
                V_imm = 0
            else:
                V_imm = 0
                z  = self.min + self.z
                z0 = z
                while z < self.z + self.max and z < z_eau:
                    V_imm += (self.section(z - z0)) * epsilon
                    z += epsilon
            Pi = V_imm * p_fluide * g
            return Pi

        
def est_une_mer(objet):
    """renvoie vrai ssi l'objet est une mer"""
    return objet.__class__.__name__ in ('mer','sea')

class film:
    
    def __init__(self, Framerate = 12, EndOfTime = 10, StartOfTime = 5):
        """crée la classe film, qui est plutôt self explanatory...
        
        monde doit être une classe world remplie d'objets(au besoin).
        FrameRate est en Hz
        EndOfTime est en seconde.
        
        Le monde est basiquement initialisé à ti = 0
        On a donc (EndOfTime * FrameRate) images de chargées"""
        global monde
        self.FR = Framerate
        self.dt = 1/ self.FR
        self.ti = StartOfTime
        self.tf = EndOfTime
        X = []
        self.T = []
        self.S = []
        self.V = []
        self.P = []
        self.Π = []
        self.objets =  []
        self.N = []
        for objet in monde.objets:
            if not est_une_mer(objet):
                X.append([])
                self.V.append([])
                self.P.append([])
                self.Π.append([])
                self.objets.append(objet)
                self.X = X
            else:
                self.N.append([]) # le niveau de la mer
        self.monde = monde
        self.done = False
    
    def realise(self):
        """réalise et mémorise le film."""
        monde = self.monde
        while monde.t < self.ti:
             monde.update(self.dt)
        self.amplitudes = []
        while monde.t < self.tf:
            monde.update(self.dt)
            P = 0
            for i in range(len(monde.objets)):
                objet = monde.objets[i]
                if not est_une_mer(objet):
                    self.X[i-P].append(objet.z)
                    self.V[i-P].append(objet.v_z)
                    self.P[i-P].append(objet.P)
                    self.Π[i-P].append(objet.Π)
                else:
                    self.N[P].append(self.monde.sea.niveau(self.objets[0].cote,self.monde.t))
                    P += 1 #compte le nombre de mers passées.
            self.T.append(monde.t)
            self.S.append(self.wave_image(monde.t))
        self.done = True
    
    def diminuer(self):
        for i in range(len([self.X,self.V,self.P,self.Π])):
            ## choix de la courbe
            courbe = [self.X,self.V,self.P,self.Π][i]
            for j in range(len(courbe)):
                data = courbe[j]
                maxi = max(data)
                mini = min(data)
                if mini == maxi:
                    mini = 0
                for k in range(len(data)):
                    valeur = data[k] - mini
                    valeur *= 1/(maxi-mini)
                    data[k] = valeur
                courbe[j] = data
            ## remise des courbes en place
            if i == 0:
                self.X = courbe
            elif i ==1:
                self.V = courbe
            elif i == 2:
                self.P = courbe
            elif i==3:
                self.Π = courbe
                
    def mer(self):
        """retourne la première mer"""
        for objet in self.monde.objets:
            if est_une_mer(objet):
                return objet
        return None
    
    
    def wave_image(self , t, cadre = [-30,30]):
        """renvoie l'image de la mer dans un cadre de dimension /cadre/ à l'instant t
        
        cadre  = [ x_min , x_max ]                                  >> Ensemble des absisses occupées par la mer simulée
        
        /!\ Ne se concentre que sur la première mer enregistré.     >> La mer construite de base dans world.__init__() avec /mer = True/
        
        /epsilon/ est l'écart spatial entre les points calculés.    >> Une valeur plus grande signifie un plus grand nombre de points, et donc un plus grand nombre de calculs."""
        y = 0
        mer = self.mer()
        S = []
        s = []
        for x in np.arange(cadre[0],cadre[-1],epsilon):
            S.append(mer.F((x,y),t))
            s.append(x)
        self.s = s
        return S
    
    def montrer(self,x,y, objet):
        """montre un cube au coordonnées désirées"""
        A = objet.cote /2 # on a un cube de coté 2A qui sera montré
        plt.plot([x+A,x+A,x-A,x-A,x+A],[y+A,y-A,y-A,y+A,y+A], color = 'k')
        plt.plot([x],[y],'r+')
        plt.text(x + 2*A, y, objet.__class__.__name__)
    
    def display(self, stat = False):
        """projete le film"""
        t = 0
        if not self.done:
            self.realise()
        for i in range(len(self.T)):
            t += self.dt
            ## écran de gauche
            plt.subplot(121,axisbg = 'w')
            plt.cla()
            for j in range(len(self.X)):
                plt.plot(self.T,self.X[j] , color = 'b')
                plt.plot(self.T,self.V[j] , color = 'g')
                plt.plot(self.T, self.P[j] , color = 'y')
                plt.plot(self.T,self.Π[j],color = 'm')
            plt.vlines(t+self.ti,0,self.max_screen(), color = 'red')
            ## écran de droite
            plt.subplot(122, axisbg = self.monde.ciel)
            plt.cla()
            for j in range(len(self.N)):
                plt.hlines(self.N[j][i],-10,+10,'g','-.')
            plt.vlines(0,-10,10,'y','-.')
            plt.plot(self.s,self.S[i])
            plt.axis('equal')
            for j in range(len(self.X)):
                self.montrer(0,self.X[j][i], self.objets[j]) # montre l'objet désiré
            plt.pause(0.1 * self.dt)
        return True
    
    def max_screen(self):
        """renvoie le point le plus haut de la courbe la plus haute du tableau de gauche"""
        maxi = 0
        for data in (self.X,self.V,self.P,self.Π):
            for  courbe in data:
                    maxi = max(courbe)
        return maxi

##algo:
if input("Voulez-vous la simulation normale?\n>>>\t") in yes:
    work = True
else:
    work = False
    
if work:
    t1 = T.clock()
    monde = world()

    if input('Voulez-vous mettre tous les objets en relations les uns aux autres?\n\t') in yes:
        Ob1 = world.cube(5, monde)
        Ob2 = world.table((5,1,1,4), monde)
        Ob3 = world.pillier((5,1,1), monde)
        Ob4 = world.pillier2((5,1,1),monde)
        Ob5 = world.pillier3((5,1,1,1,5),monde)
        Ob1.set_Mv_to(104)
        Ob2.set_Mv_to(104)
        Ob3.set_Mv_to(104)
        Ob4.set_Mv_to(104)
        Ob5.set_Mv_to(104)
    else:
        Ob  = world.cube(5, monde)
        Ob .set_Mv_to(104)
    
    F = film(Framerate = 24, EndOfTime = 15)
    F.realise()
    t2 = T.clock() - t1
    ans = None
    
    while not (ans in yes or ans in no):
        ans = input('Ready : ') 
        if ans in ('reduire','taille','r','t'):
            F.diminuer()
    t3 = T.clock() - t2
    if ans in yes:
        F.display()
        print('temps complet : ',T.clock()-t1,'\ntemps chargement : ', t2,'\ntemps attente : ',t3)
        plt.pause(2)
        plt.close()
    else :
        print(' - - - Attention, fonctionnement arrêté - - - ')
else:


    ## gen
    
    def algo_gen(n,score_lim,Surface,keep_iter = .5):
        plt.show()
        ## initialisation
        global monde
        global sample
        global score
        global keep
        global reduct
        if input("réduire?\t") in yes:
            reduct = True
        else:
            reduct = False
        if keep_iter <= 0 or keep_iter > 1:
            keep_iter = .5
        keep = keep_iter
        plt.show()
        sample = []
        monde = world()
        for i in range(n):
            sample.append(generate(Surface))
        Y1,Y2,X,x,score = [],[],[],0,[2 * score_lim]
        ## iteration
        i = n // 2
        while score[0] > score_lim:
            Y1.append( score[0] )
            Y2.append(sum(score[:i])/i)
            X.append(x)
            x += 1
            ## effets visuels
            if x%1 == 0:
                gen_show(X,Y1,Y2,sample[0])
                print(score[0],moy(score),len(sample),sep = ' | ')
            iterate = calcul(sample,Surface)
            iterate.start()
            plt.pause(2)
            iterate.join()
        ## fin
        gen_show(X,Y1,Y2,sample[0])
        plt.show()
        plt.pause(2)
        plt.close()
        return sample[0]

    class calcul(Thread):
        
        def __init__(self,sample,Surface):
            Thread.__init__(self)
            self.sample, self.Surface = sample,Surface
        def run(self):
            global sample
            global score
            sample,Surface = self.sample,self.Surface
            sample,score = iteration(sample,Surface)
    
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
    
    def objet_show(objet,e = epsilon,Couleur = 'k',alpha = 1):
        r0 = (objet.section(objet.max)/np.pi)**(.5)
        G = [0]
        D = [0]
        Z = [objet.min]
        for z in np.arange(objet.min,objet.max,e):
            r = (max(objet.S(z),e)/np.pi)**(.5)
            G.append( -r ) 
            D.append( +r )
            Z.append(  z )
        G = G + [-r0, 0 ]
        D = D + [ r0, 0 ]
        Z = Z + [objet.max,objet.max]
        if reduct:
            M = max(D)
            m = max(Z) - min(Z)
            for i in range(len(D)):
                D[i] = m*D[i]/(2*M)
                G[i] = m*G[i]/(2*M)
        plt.plot(G,Z,color = Couleur)
        plt.plot(D,Z,color = Couleur)
    
    def gen_show(X,Y1,Y2,objet):
        plt.subplot(121, axisbg = 'white')
        plt.cla()
        plt.plot(X,Y1)
        plt.plot(X,Y2)
        plt.subplot(122,axisbg = monde.ciel)
        plt.cla()
        R = (sample[:])
        for i in range(len(R)):
            objet = R[i]
            a = ((i+1) / len(R))
            try:
                objet_show(objet,Couleur = [a,a,a], alpha = (1-a)**4)
            except:
                1 + 1
            plt.pause(0.1)
        plt.pause(3)
        plt.show()
        
    
    def amplitude(objet):
        """
        Prends un objet, le situe dans le monde, puis réalise une simulation. Renvoie la valeur maxw de la vitesse /z/ enregistrée par le <film>;
        Aucun affichage graphique
        """
        global monde
        global F
        F = film(monde, Framerate = 24, EndOfTime = 15)
        F.realise()
        a = F.amplitudes[0]
        return a

    
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
        #print(f[27:])
        exec( 'global fonction_cache_124578963 \n' + f )
        return fonction_cache_124578963
        
    def iteration(sample,Surface):
        ## establishing the cost
        cost_list = []
        for f in sample:
            cost_list.append(cost(f))
        ## rangement de la chambre
        sample , cost_list = ranger(sample,cost_list)
        n = len(sample)
        i = int(n * keep)
        ## meurtre
        sample = sample [:i]
        ## reproduction
        for j in range(i):
            sample.append(modify(sample[j].clone()))
        k = 0
        while len(sample) < n:
            if np.random.random() >= 0.3:
                sample.append(modify(sample[k].clone()))
            else:
                O = modify(generate(Surface))
                for j in range( i + k ):
                    O = modify(O)
                sample.append(O)
            k += 1
            k = k % len(sample) 
        return sample,cost_list

    def child(objet1,objet2):
        """prends 2 objets et renvoie la fusion des deux (les moyennes)"""
        ## birth
        d_f1,d_f2 = objet1.function,objet2.function
        enfant    = dict()
        f1,f2,    = objet1.S, objet2.S
        for point in d_f1:
            y1 = f1(point)
            y2 = f2(point)
            enfant[point] = ( y1 + y2 ) / 2
        for point in d_f2:
            y1,y2 = f1(point),f2(point)
            enfant[point] = (y1 + y2) / 2
        ## rangement
        return world.defineable(monde,enfant)

    def moy(iterable):
        """returns an avarged number out of the iterable"""
        sum = 0
        for truc in iterable:
            sum += truc
        return sum / len(iterable)
    
    def ranger(liste, valeur):
        """renvoie le couple/liste/,/valeur/ rangés pour valeur croissante"""
        for i in range(len(liste)):
            item = {"value" : valeur [i], "data" : liste[i] }
            for j in range(i,len(liste)):
                if valeur[j] < item["value"]:
                    item["value"],valeur[j] = valeur[j],item["value"]
                    item["data" ],liste [j] = liste [j],item["data" ]
            valeur[i] = item["value"]
            liste [i] = item["data" ]
        return liste,valeur
    
    def generate(Surface):
        """ renvoie un objet légal """
        global monde
        f = creer(a = 0 ,b = 20,f_a = 0,f_b = Surface)
        Ob = world.defineable(monde,f,0,20)
        return Ob
    
    def cost(objet):
        """Fonction C qui à chaque fonction 
        R  -> R
        f: z |-> Section(z)
        renvoie une valeur basée sur l'amplitude du mouvement subit par l'objet définit par f et un monde sans argument"""
        a,b = v_max(objet) , amplitude(objet)
        return a**2 + b**2
    
    def modify(objet):
        """renvoie une version légèrement modifiée de l'objet"""
        if np.random.randint(20) > 5 :
            f = objet.function
            f = warp(f)
            objet.mutate(f)
        else:
            try:
                return child(objet,sample[0])
            except:
                return child(objet,generate(objet.function[objet.max]))
        
        return objet
    
    def warp(f):
        """prends le dictionnaire et sois modifie la valeur, soit en ajoute une autre.
        
        Ne modifie jamais la seconde valeur, qui sera utilisée afin de stocker l'amplitude"""
        r = np.random.randint(1,3)
        F = fonction(f)
        if len(f) <= 2:
            r = 1
        if r == 1:
            ##make new points
            i = 0
            A = min(f)
            B = max(f)
            R = np.random.random()
            R = A + (B-A) * R
            r = np.random.random()
            rand = (F(R) + r * np.random.normal())/r
            f[R] = rand
        else:
            ## move old ones
            i = np.random.randint(0,len(f))
            j = 0
            dont = False
            for truc in f:
                if j == i :
                    if truc != max(f):
                        A = truc
                    else:
                        dont = True
                j += 1
            if not dont:
                a = f[A]
                v = (abs(a) + 1) ** 0.47                    # on fait plus varier les valeurs plus grandes, d'où une f croissante
                f[a] = a + (np.random.random() * 2 - 1) * v
        return f
    
    def amplitude(objet):
        return objet.amplitude(ti = 2,tf = 16)
    
    def v_max(objet):
        return objet.v_max(ti = 2,tf = 16)