#!/usr/bin/env ruby
# encoding: binary
#
# HIST
# - 06/12/2022: 
#   - amélioration génération du proto pour les élèves
#   - création d'une méthode do_gitpush (qui envoie le projet vers tous les 
#     serveurs distants déclarés dans ce projet git)
#     Effectue par défaut la commande suivante pour tous les remotes
#        git push --all --tags myremote
#   - ajout d'une target git_push (qui envoie le projet vers tous les 
# - 01/12/2021: ajout cible relatives à la génération de la doc
# - 23/02/2020: Amélioration de la cible distclean(alias dc)
#     pour pouvoir supprimer 11000 fichiers du répertoire "_tmp/"
#     La commande rm (ou le shell ?) ne peut pas supprimer 11000 fichiers
#     du répertoire  "_tmp/" !`
#     Je passe donc par ruby pour supprimer les fichiers par groupe de 1000
# - xx/xx/201x: nombreuses modifications
# - xx/xx/2014: création pour corriger le projet de Alain Faye
#

require 'colorize'
require 'amazing_print'

# Retourne le répertoire contenant ce Rakefile (quelque soit l'endroit où l'on est)
# e.g.  /Users/diam/live/public_html/uma/pub/work/sf/mpro_proj/Rakefile
def appdir
  return "#{File.dirname(__FILE__)}"
end

PROJET = "seqata"

# Constante spécifique personnelle (non portable)
TOPDIR = File.dirname(File.dirname(appdir))

# Je choisis de générer la même doc pour le projet `Seqata.jl` ou le `proto_seqata`
DESTDOC = "#{TOPDIR}/www_minisme_seqata/docs/seqata_docs/"


task :default => [:help]


# Crée une archive "-gitarc.txz" datée du répertoire dir à coté de dir
# dir: répertoire de l'archive à créer
# suf: suffixe à ajouter dans le nom de l'archive datée (apres la date
#      et avant le suffixe réel
#
# Vide le répertoire
#
# Voir de mon script shell ~/local/bin/z (pour modèle) :
#
#   bdir=$(basename $(pwd))
#   new_bdir=$bdir-`dateString`
#   git archive --prefix "${new_bdir}-gitarc/" ${what} \
#       | xz -c  > ../${new_bdir}-gitarc.txz
#
def do_gitarc(dir, suf="")

  if !Dir.exist?("#{dir}/.git")
    puts "Ce projet n'est pas géré par git"
    puts "   dir=#{dir}"
    puts "Opération annulée."
    exit(0)
  end

  # Paramètre spécifique à git
  what = "HEAD:"

  # On s'assure que si le suffixe existe, il commence par "-"
  suf = "-"+suf unless suf.start_with?("-")

  # On se positionne à la racine du projet (même si pas nécessaire avec git)
  dir_ori = pwd
  cd appdir

  bname = File.basename(dir)
  date = Time.now.strftime("%Y%m%d_%Hh%M")
  datename = "#{bname}-#{date}-gitarc"

  cmd  = ""
  cmd += " git archive --prefix \"#{datename}/\"  #{what}"
  cmd += " | xz -c  > \"../#{datename}#{suf}.txz\" "

  # puts "POUR INFO : cmd=#{cmd}"
  sh cmd
  # puts "=> fait: #{cmd}"

  cd dir_ori
end

# Génère la doc automatique
# - génère la doc par make.jl
# - transfert le build dans le sous-répertoire distrib
# - met à jour le lien vers la dernìere version
# - Supprime les anciennes doc de distrib pour ne conserver que la dernìère
def do_docs()
  do_docs_make
  do_docs_clean_build
end

# Enchaine Génèration, nettoyage et déploiement de la doc sous le public_html
def do_docs_push()
  do_docs_push_to_local
  do_docs_push_to_minisme
end

# Enchaine Génèration, nettoyage et déploiement de la doc sous le public_html
def do_docs_prod()
  starttime = Time.now
  do_docs
  do_docs_push
  puts "Durée de la génération de la doc : #{(Time.now - starttime).round(2)}s"
end


# Génère la doc automatique
# - génère la doc par make.jl
# - transfert le build dans le sous-répertoire distrib
# - met à jour le lien vers la dernìere version
# - vide le répertoire docs/build
#
def do_docs_make()

  puts "\nGénère la doc datée et met à jour le lien vers la dernière version"

  t1 = Time.now

  date = Time.now.strftime("%Y%m%d_%Hh%M")

  # On veut récupérer (e.g) la chaine "seqata" à partir du nom du répertoire `Seqata.jl`
  bname = File.basename(appdir, ".jl").downcase

  # Création de l'auto-doc par julia
  # sh "julia --color=yes   #{appdir}/docs/make.jl"
  sh "#{appdir}/docs/make.jl"
  # recopie du .htaccess (pour redirection vers les fichiers index.html)
  sh "cp -a  \"#{appdir}/docs/.htaccess_allow_access_but_listing_http24\" \
             \"#{appdir}/docs/build/.htaccess\""

  # récopie de la doc générée vers un répertoire daté pour la distrib
  sh "cp -a  \"#{appdir}/docs/build\" \"#{appdir}/distrib/#{bname}_docs-#{date}\""

  # on recrée un lien sans date vers la dernìere version de cette doc
  sh "rm -f \"#{appdir}/distrib/#{bname}_docs\""
  sh "cd \"#{appdir}/distrib/\" &&  ln -s \"#{bname}_docs-#{date}\" \"#{bname}_docs\""
  # sh "ln -s \"#{appdir}/distrib/#{bname}_docs-#{date}\" \"#{appdir}/distrib/#{bname}_docs\""

  puts "\nDocs générée en #{(Time.now-t1).round(1)}sec : \"#{bname}_docs-#{date}\" "

  puts "\nSupprime les anciennes doc de distrib pour ne conserver que la dernière"
  do_clean_old_docs
end


# Supprime les anciennes doc de distrib pour ne conserver que la dernìère
def do_clean_old_docs

  puts "\nSupprime les anciennes doc de distrib pour ne conserver que la dernìère"

  bname = File.basename(appdir, ".jl").downcase
  gpat =  "#{bname}_docs-202?????_??h??"
  dirs =  Dir.glob("#{appdir}/distrib/#{gpat}")
  dirs.sort! # En principe fait par défaut par Dir.glob

  lastdir = dirs.pop
  puts "on ne conserve que : \"#{lastdir}\""

  for dir in dirs
    sh "rm -rf \"#{dir}\""
  end
end

# Supprime les anciennes doc de distrib pour ne conserver que la dernìère
# et supprime le comteni du répertoire docs/build
#
def do_docs_clean_build
    puts "\nSupprime le contenu de docs/build/"
    # sh "rm -rf \"#{appdir}/docs/build/*\""  # le "*" serait masqué pour le shell
    sh "rm -rf \"#{appdir}/docs/build/\"" + "*"
    sh "rm -rf \"#{appdir}/docs/build/\"" + ".htaccess"
end

def do_tmpclean
    puts "\nSupprime le contenu de _tmp/ (par paquets de 1000 fichiers)"
    dir="_tmp"
    # La commande rm (ou le shell ?) ne peut pas supprimer 11000 fichiers !
    # Je passe donc par ruby pour supprimer les fichiers par groupe de 100`
    #
    # sh "rm -f #{dir}/*"   # Ceci fonctionnait avec un nombre raisonnable du fichiers

    Dir["_tmp/*"].each_slice(1000) do |files|
        puts "files.length=#{files.length}"
        # puts "files[0]=#{files[0]}"
        # sh "ls -altr #{files.join(" ")}", verbose: false
        sh "\\rm -- #{files.join(" ")}", verbose: false
    end
end

def do_distclean
    do_tmpclean
    do_docs_clean_build

    puts "\nSupprime le contenu de distrib"
    sh "rm -rf \"#{appdir}/distrib/\"" + "*"
end

# Copie la sous-rép distrib dans le sous-répertoire web local
# Le répertoire dst **doit** exister
def do_docs_push_to_local

  if ENV['USER'] != "diam"
    puts "\nDésolé, cette action n'est utilisable que par l'utilisateur \"diam\""
    exit
  end

  src = "#{appdir}/distrib/"
  dst = DESTDOC + "/"

  puts "\nCopie le sous répertoire distrib vers #{dst}"
  # sh "cp -n -a  \"#{src}\"  \"#{dst}\""
  sh "rsync -av --stats --delete \"#{src}\"  \"#{dst}\""

end

# Envoie le sous-répertoire web local vers le serveur minisme
# (en vue de son accès depuis le site minisme)
def do_docs_push_to_minisme

  if ENV['USER'] != "diam"
    puts "\nDésolé, cette action n'est utilsable que par l'utilisateur \"diam\""
    exit
  end

  sh "minisme push_seqata"

  # puts "\nEnvoye du répertoire web local vers le serveur".red
  # puts "vous devez saisir MANUELLEMENT l'alias suivant \"Seqata_push_to_minisme\"".red
  # puts "   Seqata_push_to_minisme".red
  # sh "Seqata_push_to_minisme"
end


# Crée un archive "txz" datée du répertoire dir à coté de dir
# dir: répertoire de l'archive à créer
# suf: sufixe à ajouter dans le nom de l'archive datée (apràs la date
#      et avant le suffixe réel
def do_txz(dir, suf="")
  dir_ori = pwd
  dir_parent = File.dirname(dir)
  bname = File.basename(dir)
  date = Time.now.strftime("%Y%m%d_%Hh%M")
  ## datename = "#{bname}—#{date}" ##!!!  BUG du "-" !!!
  datename = "#{bname}-#{date}"

  cmd = "cd #{dir_parent}"
  cmd += " && cp -Rp '#{bname}' '#{datename}' "
  cmd += " && tar cf - #{datename}"
  cmd += " | xz > '#{datename}#{suf}.txz' &&  rm -Rf '#{datename}'"
  cmd += " && cd #{dir_ori}"

  # puts "POUR INFO : cmd=#{cmd}"
  sh cmd
  # puts "=> fait: #{cmd}"
end

# Crée un archive "tbz" datée du répertoire dir à coté de dir
# dir: répertoire de l'archive à créer
# suf: sufixe à ajouter dans le nom de l'archive datée (apràs la date
#      et avant le suffixe réel
def do_tbz(dir, suf="")
  dir_ori = pwd
  dir_parent = File.dirname(dir)
  bname = File.basename(dir)
  date = Time.now.strftime("%Y%m%d_%Hh%M")
  ## datename = "#{bname}—#{date}" ##!!!  BUG du "-" !!!
  datename = "#{bname}-#{date}"

  cmd = "cd #{dir_parent}"
  cmd += " && cp -Rp '#{bname}' '#{datename}' "
  cmd += " && tar cf - #{datename}"
  cmd += " | bzip2 > '#{datename}#{suf}.tbz' &&  rm -Rf '#{datename}'"
  cmd += " && cd #{dir_ori}"

  # puts "POUR INFO : cmd=#{cmd}"
  sh cmd
  # puts "=> fait: #{cmd}"
end

# Crée un archive "zip" datée du répertoire dir à coté de dir
# dir: répertoire de l'archive à créer
# suf: sufixe à ajouter dans le nom de l'archive datée (apràs la date
#      et avant le suffixe réel
def do_zip(dir, suf="")
  dir_ori = pwd
  dir_parent = File.dirname(dir)
  bname = File.basename(dir)
  date = Time.now.strftime("%Y%m%d_%Hh%M")
  ## datename = "#{bname}—#{date}" ##!!!  BUG du "-" !!!
  datename = "#{bname}-#{date}"

  cmd = "cd #{dir_parent}"
  cmd += " && cp -Rp '#{bname}' '#{datename}' "
  cmd += " && zip -r -y -o -q -9 '#{datename}#{suf}.zip'  '#{datename}'"
  cmd += " && rm -Rf '#{datename}'"
  cmd += " && cd #{dir_ori}"

  # puts "POUR INFO : cmd=#{cmd}"
  sh cmd
  # puts "=> fait: #{cmd}"
end

# Crée une archive datée et réduite d'une branche `proto` pour les élèves
# - branch: la branche git à archiver qui doit-être active et clean
#           (e.g proto_p2022).
# - arc_bname: le nom de base de l'archive (e.g proto_seqata_p2022)
#
# Fait un git push vers l'ensemble des serveurs distants
# - gère l'options booléenne all (pour ajouter --all à git push)
# - gère l'options booléenne tags (pour ajouter --tags à git push)
# - gère une option dry (défaut true)
# 
# TODO
# - voir pour gérer l'option -u dans "git push -u" (quelle valeur par défaut)
#
def do_gitpush(all:true, tags:true, dry:false)
  #
  # 1. on vérifie que lq branche courante est clean (i.e. commitée)
  #
  status = %x(git status)
  if !status.match? /nothing to commit, working tree clean/
    puts "Erreur : le branche git active n'est pas clean !".red
    exit(1)
  end
  puts "Le branche git courante est bien clean".green

  #
  # 2. extraction des remote serveurs
  #
  remotes = %x(git remote).scan(/\S+/)

  #
  # 3. construcion de la liste des commandes
  #
  cmds = []
  for rem in remotes
    cmds << "git push --all #{rem}"     if all
    cmds << "git push --tags #{rem}"    if tags
    cmds << "git push #{rem}"           if !all && !tags
  end

  #
  # 4. construcion de la liste des commandes
  #
  for cmd in cmds
    if dry
      puts "DRY: #{cmd}"
    else
      # puts "EXCECUTION (dry=false): #{cmd}"
      sh cmd
      puts "#{cmd} FAIT".green
    end
  end
end


# Crée une archive datée et réduite d'une branche `proto` pour les élèves
# - branch: la branche git à archiver qui doit-être active et clean
#           (e.g proto_p2022).
# - arc_bname: le nom de base de l'archive (e.g proto_seqata_p2022)
#
def do_proto(branch:"proto_p2023", arc_bname:"proto_seqata_p2023")
  if ENV['USER'] != "diam"
    puts "\nDésolé, cette action n'est utilisable que par l'utilisateur \"diam\""
    exit
  end

  #
  # 1. Vérification de la branche à extraire (active et clean ?)
  #
  branches = %x(git branch)
  m = branches.match(/\* ([-\w]+)/)
  # ap m
  if !m
    puts "Erreur : impossible de trouver la branche git active".red
    exit(1)
  end
  if m[1] != branch
    puts "Erreur : le branche git active \"#{m[1]}\" devrait être \"#{branch}\" !".red
    exit(1)
  end
  puts "On est bien sur la branche git \"#{branch}\"".green

  status = %x(git status)
  if !status.match? /nothing to commit, working tree clean/
    puts "Erreur : le branche git active \"#{branch}\" n'est pas clean !".red
    exit(1)
  end
  puts "Le branche git \"#{branch}\" est bien clean".green

  puts "Nettoyage des fichiers temporaires (par rake disclean)...".yellow
  do_distclean
  puts "Nettoyage des fichiers temporaires (par rake disclean)... FAIT".green

  #
  # 2. Création d'une copie datée dans le répertoire parent
  #

  date = Time.now.strftime("%Y%m%d_%Hh%M")
  arc_dated_name = "#{arc_bname}-#{date}"
  srcdir = appdir
  dstdir = File.dirname(appdir) + "/" + arc_dated_name
  puts "Création du repertoire \"#{dstdir}\"... ".yellow
  cmd = "cp -Rp '#{srcdir}' '#{dstdir}'"
  sh cmd
  puts "Création du répertoire \"#{dstdir}\"... FAIT".green

  #
  # 3. Suppression de sous-répertoire ne concernant pas les élèves
  #

  # Liste de paternes à supprimer.
  # Attention à l'utilisation de la commande `rake sh`` :
  # Pour utiliser les métacaractères dans la commande rake `sh`, il faudra lui
  # passer une String et non pas un Array
  # - On ne peut pas faire :   sh ["ls", "*.md"]
  # - Il faudra faire plutot :        sh "ls *.md"
  # Conclusion : éviter les métacaractères quand c'est possible
  pats_to_delete = %w(
    ADMIN
    ADMIN_SEQATA
    .git
    sols/alp_{0[2-9],1[0-3]}_p*.sol

    src/dynprog_timing_solver.jl
    src/faye_timing_solver.jl

    src/taboo_solver_FROM_ALAP.jl
    src/grasp_solver.jl
    src/greedy_solver.jl
    src/main_ant.jl
    src/main_grasp.jl
    src/main_greedy.jl
  )
  puts "Suppression de certains répertoires avant compression...".yellow
  cd dstdir # DANGER: NE PAS OUBLIER CETTE COMMANDE !!
  for pat in pats_to_delete
    # cmd = "echo ls -al #{pat}"
    cmd = "rm -rf #{pat}" # DANGER!
    puts "todo cmd: #{cmd}".yellow
    sh cmd
  end
  puts "Suppression de certains répertoires avant compression... FAIT".green

  #
  # 4. Compression de l'archive
  #
  cd File.dirname(dstdir)
  puts "Création de l'archive compressée...".yellow
  sh "zip -r -y -o -q -9 '#{arc_dated_name}.zip'  '#{arc_dated_name}'"
  puts "Création de l'archive compressée... FAIT".green

  #
  # 5. Suppression de la copie intermédiaire
  #
  puts "Suppression du répertoire temporaire...".yellow
  cmd = %W(rm -r #{arc_dated_name}) # ni métacaractère ni espace dans arc_dated_name
  # ap cmd
  sh *cmd
  puts "Suppression du répertoire temporaire... FAIT".green

  cd appdir  # pour revenir au répertoire de départ à la fin

  puts "do_proto : fait !".green
end


#
# Quelques cibles standard (compression...)
#

desc "Crée une archive datée (xxx.txz) dans le répertoire parent"
task :txz do
  do_txz pwd
end
desc "Crée une archive datée (xxx-stamp.txz) dans le répertoire parent"
task :txzstamp do
  do_txz pwd, "-stamp"
end
desc "Crée une archive datée (xxx.tbz) dans le répertoire parent"
task :tbz do
  do_tbz pwd
end
desc "Crée une archive datée (xxx-stamp.tbz) dans le répertoire parent"
task :tbzstamp do
  do_tbz pwd, "-stamp"
end
desc "Crée une archive datée (xxx.zip) dans le répertoire parent"
task :zip do
  do_zip pwd
end
desc "Crée une archive datée (xxx-stamp.zip) dans le répertoire parent"
task :zipstamp do
  do_zip pwd, "-stamp"
end

#
# Quelques cibles spécifiques au projet courant
#


desc "Supprime **tous les fichiers** du sous-répertoire générées _tmp/"
task :tmpclean do |t|
  do_tmpclean
end

desc "Supprime **tous les fichiers** du sous-répertoire générées _tmp/"
task :distclean do |t|
  do_distclean
end
# task :dc, [:par1, :par2] => :distclean
task :dc => :distclean

desc "Crée une archive datée (xxx-gitarc.txz) avec git " +
      'eg : rake gitarc suf="avec_descent" ' +
      "(n'est disponible que pour un projet sous git)"
task :gitarc do
  suf = ENV["suf"] || ""
  do_gitarc pwd, suf
end


desc "Génère la doc automatique dans distrib et nettoie anciennes version de la doc"
task :docs do
  do_docs
end

desc "déploie le sous-rép distrib vers les répertoires web local et distant (diam only)"
task :docs_push do
  do_docs_push
end


desc "Déploie le répertoire distrib vers le web local et vers le serveur (diam only)"
task :docs_prod do
  do_docs_prod
end

desc "Crée une archive datée et réduite de la branche proto (diam only)"
task :proto do
  do_proto
end

desc "EN COURS DÉVELOPPEMENT: push vers l'ensemble des serveurs distant"
task :gitpush do
  # do_gitpush dry:true, all:false, tags:false   # DRY POUR TEST
  # do_gitpush dry:true, all:true, tags:true     # DRY POUR TEST
  # do_gitpush dry:false, all:true, tags:true    # valeur par défaut
  do_gitpush
end

desc "Fournit une aide minimaliste de ce Rakefile"
task :help do
  puts ""
  puts "Rakefile spécifique au projet \"#{PROJET}\""
  puts ""
  puts "Quelques options utiles de rake "
  puts ""
  puts "   rake -T  : liste des taches documentées par desc "
  puts "   rake -P  : liste des dépendances"
  puts "   rake -D  : Describe"
  puts ""
  puts "Quelques exemples"
  puts %Q(
    rake docs         # génère la doc automatique datée dans distrib/xxx  (avec lien sym)
    rake docs_push    # pousse le contenu de distrib vers site web locale et distant (diam only)
    rake docs_prod    # génère puis pousse la doc vers serveur web
    rake tmpclean     # supprime tous les fichiers générés dans _tmp/
    rake dc           # alias pour disclean
    rake distclean    # supprime tout répertoire rempli automatiquement
                      # (_tmp, distrib, docs/build)
    rake zip ou tbz ou txz # crée une archive datée dans répertoire parent

    rake gitarc       # archive git datée et avec label (exemple : rake gitarc suf="prefinale")
    rake gitpush      # effectue un "git push" sur tous les serveurs remote déclarés
    rake -T           # affichées toutes les cibles disponibles
  ).gsub /^ {4}/, ""
end
