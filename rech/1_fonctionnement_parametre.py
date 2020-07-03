# coding: utf-8
 
from tkinter import *
import webbrowser 
import tkinter.filedialog
from tkinter import ttk

  
def aide_cours():
    webbrowser.open_new("http://www.patricebacchin.fr/cours/membrane/")

def enregistrer_fichier(mesFormats,interface):
    tkinter.filedialog.asksaveasfile(title="Enregistrer sous … un fichier",filetypes=[('CSV files','.csv')])
    
def ouvrir_fichier():
    tkinter.filedialog.askopenfilename(title="Ouvrir un fichier",filetypes=[('CSV files','.csv')]) 

def sauv_donnée():
    Fichier = open("Valeur rentrée.txt", 'w')
    Fichier.writelines(str(L))
    Fichier.close()
    
def liste():
    L=[]
    global const, arg4
    var=arg4.get()
    print(var)
    L.append(var)
    print(L)

def paramètre(arg1, arg2, arg3):
    global const, arg4
    arg4 = DoubleVar()
    text = Label(arg1, text=arg2, width=10)
    entree=Entry(arg1, textvariable=arg4, width=30)
    text.grid(row=arg3, column=0)
    bouton=Button(arg1, text="Entrée", command=liste)
    entree.grid(row=arg3, column=1)
    bouton.grid(row=arg3, column=2)
    const += 1

def onglet(arg1, arg2):
    n = ttk.Notebook(arg1)
    n.grid(row=5, column=0)
    o1 = Frame(n)
    o1.grid(row=5, column=0)
    n.add(o1, text = arg2)  
    
def données():
    fen_1 = Frame(fen)
    fen_1.grid(row=5, column=0, columnspan=3, sticky="nsew", pady=10)
    onglet(fen_1, "Calcul nombre de Reynolds")
    paramètre(fen_1, "Densité (noté {0:^1s})".format(chr(961)), 6)
#    paramètre(fen_1, "Vitesse moyenne de la fibre", 7, "val2")
#    paramètre(fen_1, "Viscosité dynamique (noté {0:^1s})".format(chr(956)), 8, "val3")
#    paramètre(fen_1, "Rayon du papier", 9, "val4")
    button_1 = Button(fen_1, text="Valider", font=("Times New Roman", 25), bg='white', fg='#C1B5B3', command=sauv_donnée)
    button_1.grid(row=6, column=3, pady=10)

#création d'une première fenètre
fen = Tk()

#personnalisation de cette fenêtre
fen.title("Simulation du procédé membranaire")
fen.geometry("1080x720")
fen.minsize(480,360)
fen.iconbitmap("")
fen.config(background='#C1B5B3')

# Configuration du gestionnaire de grille
fen.rowconfigure(0, weight=3)
fen.columnconfigure(0, weight=3)

#Création de la barre d'outil
menubar = Menu(fen)
fen.config(menu=menubar)
menufichier = Menu(menubar,tearoff=0)
menubar.add_cascade(label="Fichier", menu=menufichier)
menufichier.add_command(label="Ouvrir  ", command=ouvrir_fichier)
menufichier.add_separator()
menufichier.add_command(label="Enregistrer")
menufichier.add_command(label="Enregistrer sous", command=enregistrer_fichier)
menufichier.add_separator()
menufichier.add_command(label="Quitter", command=fen.destroy) 

#Programme principal
const = 0
arg4 = 0
#Créer la frame
frame =Frame(fen, bg='#C1B5B3')

#Ajout texte pour rentrer à la page d'acceuil
label_title = Label(frame, text="Accueil", font=("Times New Roman", 40), bg='#C1B5B3', fg='white')
label_title.grid(row=1, column=0, columnspan=2,sticky="nsew", pady=10)

#Ajouter un premier bouton
button_1 = Button(frame, text="Données utiles pour la simulation", font=("Times New Roman", 25), bg='white', fg='#C1B5B3', command=données)
button_1.grid(row=2, column=0, columnspan=2,sticky="nsew", pady=10)

#Ajouter un troisième bouton
button_1 = Button(frame, text="Commencer la simulation", font=("Times New Roman", 25), bg='white', fg='#C1B5B3')
button_1.grid(row=3, column=0, columnspan=2,sticky="nsew", pady=10)

#Ajouter un second bouton
button_3 = Button(frame, text="Aide", font=("Times New Roman", 25), bg='white', fg='#C1B5B3', command=aide_cours)
button_3.grid(row=4, column=0, columnspan=2,sticky="nsew", pady=10)

#Ajouter
frame.grid(row=0, column=0, columnspan=2)

#Affichage de la fenetre
fen.mainloop()
