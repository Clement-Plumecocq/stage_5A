# Template_Rapport_Stage_Polytech_Latex
Template de rapport de stage (ou autre) en LateX utilisable à Polytech Marseille ou dans n'importe quelle école du réseau Polytech sous réserve de changer les images.
Nom du PDF :

__AAAA_MM_5A_SPECIALITE_NOM_PRENOM__

Comment l'utiliser ?
Télécharger l'archive ou cloner le repo :

<ins>Remarque :</ins> Je recommande ici d'installer TexLive. En effet, les logiciels Texstudio, VS Code et Sublime-Text-3 sont des éditeurs de texte, mais votre PC a besoin d'un compilateur LaTeX pour fonctionner.

# Configuration TexStudio
- Dans **Compilation**, choisir le paramètre de compilation PDFLaTeX : ``pdflatex.exe -synctex=1 -interaction=nonstopmode -shell-escape %.tex``
- Dans **Production**, choisir Biber en tant que moteur de bibliographie par défaut.
- Bien penser à appuyer sur F9 pour mettre à jour le glossaire quand il est modifié avant de compiler.
- Penser également à créer un user command "Make Nomenclature" : ``makeindex -s nomencl.ist -t %.nlg -o %.nls %.nlo`` (https://www.youtube.com/watch?v=kW97Yv0-QC4)
- Compiler avec F5 pour compiler et visualiser. La compilation peut se lancer plusieurs fois pour prendre en compte la bibliographie, bien attendre qu'il soit indiqué : ``system returned with code 1 Processus terminé normalement``. Autrement, vérifier les logs et corriger les erreurs.
