#!/usr/bin/env ruby

# PRÉREQUIS POUR INSTALLER CETTE COMMANDE
# - installer ruby
# - ajourer les gems (i.e packages) suivants
# gem install amazing_print colorize json optimist logger ostruct fileutils abbrev
#
# Ce code permet d'effectuer des tests fonctionnels pour  des applications en
# Ligne de Commande (CLI pour Command Line Interface)
#
# Les tests sont décrits dans un fichier JSON dont le nom par défaut est
# ./Funtests.json ou ./test/Funtests.json
#
# Exemple :
#
# [
# {
#     "sid": "help_h",
#     "desc": "Version minimaliste (test de la gestion des arguments)",
#     "cmd": "./bin/run.jl -h",
#     "res_desc": "Affiche l'aide en ligne du programme",
#     "exp_rpat": "(?m:usage:.*--outdir OUTDIR.*--version)",
#     "exp_duration": 9,
#     "exp_duration_margin": 3,
#     "exp_length": 8462,
#     "exp_length_margin": -1
# },
#   ...
# {
#     ...
#
# }]
#
# TODO
# - prévoir affichage du temps total prévu (somme des exp_durations)
#
#
# Dépendance à ajouter sous unix :
#    gem install awesome_print colorize json optimist logger ostruct
#
# require 'awesome_print'
require 'amazing_print'
require 'colorize'
require 'json'
require 'optimist'
require 'logger' # utilisé en interne par Madi.run_command
require 'ostruct'
require 'fileutils'
require 'abbrev'

# M'est plus utilisé depuis le 26/05/2020 (remplacé par --level => args.level)
# $log = Logger.new(STDOUT)
# $log.level = :warn # PAS ICI !

$EXE_NAME = File.basename(__FILE__)

# Quelques commandes perso pour tester et valider une nouvelle distribution de Julia et des ses packages.

# Décrit un test complet avec les attributs suivants (À COMPLÉTER):
#
# sid: identifiant si possible unique du test (chaine de caractère)
# cmd: ommande unix à exécuter pour le test
# desc: ligne de description du test
# res_desc:  description détaillée du résultat
# exp_rpat: pattern regexp. peut-être "" (pas de "//"
#    Dans le json, es option de regexp pouvent être passer à l'intérieur de la
#    patterne
#   "exp_rpat": "(?i:VALIDATEUR/SOLVER POUR (ALAP|SEQATA))",
#   "exp_rpat": "(?m:usage:.*--outdir OUTDIR.*--version)",
#
# exp_duration : durée attendu de l'exécution en secondes (-1 pour ignorer)
# exp_duration_margin : marge d'erreur sur la durée en secondes (-1 pour infini)
# exp_length : londueur de la chaine de sortie en octets (-1 pour ignorer)
# exp_length_margin : marge d'erreur sur la chaine de sortie (-1 pour infini)
#
# Le constucteur accepte :
# - soit un tableau des argument
# - soit une chaine au format JSON
#
class Test
  attr_accessor :sid, :desc, :res_desc, :cmd, :exp_rpat
  attr_accessor :exp_duration, :exp_duration_margin, :exp_length, :exp_length_margin
  @sid = "" # string id (si possible unique)
  @desc = ""
  @cmd = ""
  @res_desc = ""
  @exp_rpat = "" # regex de validation du résultat
  @exp_duration = -1 # secondes  (-1 pour non testé)
  @exp_duration_margin = -1 # secondes ; -1 pour infini
  @exp_length = -1 # octets (-1 pour non testé)
  @exp_length_margin = -1 # octets ; -1 pour infini

  def initialize array_or_json
    if array_or_json.is_a? String
      specs = JSON.parse(array_or_json)
    elsif array_or_json.is_a? Hash
      specs = array_or_json
    else
      raise "Test type du constructeur de Test inconnu : typof(array_or_json)"
    end
    @sid = specs["sid"]
    @desc = specs["desc"]
    @cmd = specs["cmd"]
    @res_desc = specs["res_desc"] if specs["res_desc"]
    @exp_rpat = Regexp.new(specs["exp_rpat"]) if !specs["exp_rpat"].empty?
    @exp_duration = specs["exp_duration"] if specs["exp_duration"]
    @exp_duration_margin = specs["exp_duration_margin"] if specs["exp_duration_margin"]
    @exp_length = specs["exp_length"] if specs["exp_length"]
    @exp_length_margin = specs["exp_length_margin"] if specs["exp_length_margin"]
  end

  def to_s
    # return "@#{@sid}: #{@desc}".blue.bold.on_light_yellow
    return "@#{@sid}: #{@desc}".blue.underline
  end

  def to_s_long(dump: false)
    # txt = "@#{@sid}: #{@desc}".blue.bold.on_light_yellow << "\n"
    txt = "@#{@sid}: #{@desc}".blue.underline << "\n"
    txt << "  cmd:        #{@cmd}\n"
    txt << "  res_desc:   #{@res_desc}\n"       if dump || !@res_desc.empty?
    txt << "  exp_rpat:   #{@exp_rpat}\n"       if dump || @exp_rpat.class == Regexp
    txt << "  exp_duration: #{@exp_duration}\n" if dump || @exp_duration != 0
    txt << "  exp_duration_margin: #{@exp_duration_margin}\n" if dump || @exp_duration_margin != -1
    txt << "  exp_length: #{@exp_length}\n"     if dump || @exp_length != 0
    txt << "  exp_length_margin: #{@exp_length_margin}\n"     if dump || @exp_length_margin != -1
    return txt
  end
end

# get_funtests: retourne le Vector des objets Test à traiter
#
# funtestsfile : fichier json contenant les test (e.g Funtests.json)
# testpats : tableau de regex pour sélectionner un sous-ensemble de tests
#
def get_funtests(funtestsfile, testpats: [] )

  json_txt = File.read(funtestsfile)
  specs = JSON.parse(json_txt)

  if  testpats.empty?
    testpats = [".*"]
  end

  funtests = Array.new
  for spec in specs
    # On ignore les entrées qui ne disposent pas d'un champs cmd
    # ce qui permet d'ajouter de pseudo commentaire dans le json
    next unless (spec["sid"] && spec["cmd"])
    funtests.push(Test.new(spec))
  end

  # On filtre par rapport au champ sid
  subfuntests = funtests.select do |test|
    for pat in testpats
      res = false
      if test.sid.match?(pat)
        res = true
        break
      end
    end
    res
  end

  # Si on a rien trouvé, on filtre par rapport au champs desc
  if subfuntests.empty?
    subfuntests = funtests.select do |test|
      for pat in testpats
        res = false
        if test.desc.match?(pat)
          res = true
          break
        end
      end
      res
    end
  end
  subfuntests
end

# CRITIQUES DU GEM Optimist
# - les options semblent immutables, ce qui oblige à passer par un attribut
#   supplémentaire
#   Exemple  @gopts.funtestsfile non modifiable, Je crée donc @funtestsfile
class Args
  attr_reader :gopts, :command, :copts
  attr_reader :testpats, :funtestsfile, :level
  @@commands = {
    'help'  => "Affiche aide et exemples",
    'list'  => 'liste les tests disponibles',
    'proto' => 'génère un fichier d\'exemple "Funtests.json"',
    'run'   => 'exécute les tests (tous ou ceux en argument)',
  }

  def Args.show_commands
    puts Args.to_s_commands
  end

  def Args.to_s_commands
    txt = ""
    @@commands.each do |cmd, desc|
      txt << "  %-10s ".blue.bold % [cmd]
      txt << desc << "\n"
    end
    txt
  end

  def initialize
    @testpats = []
    @funtestsfile = "UNDEF"
    @level = 1
    @copts, @command = nil, nil
    @gopts = Optimist::options do
      version "#{$EXE_NAME} v0.23-20221117"
      banner "Usage:"
      banner "  #{$EXE_NAME} [options] [<command> [suboptions]]\n "
      banner "Listes des sous-commandes :"
      banner Args.to_s_commands
      banner ""
      banner "Options:"
      opt :version, "Affiche la version", :short => "-V"
      opt :help, "Affiche la syntaxe", :short => "-h"
      opt :show_opts, "Affiche la valeur actuelle des options", :short => "-p"
      opt :funtestsfile, "fichier json des tests (défaut: Funtests.json)",
          :short => "-f", :default => "UNDEF"
      stop_on_unknown
    end

    @command = ARGV.shift
    @command = unabbrev(@command)
    if !@command
      # manque la sous-commande.
      # Soit on impose une commande par défaut par :
      # @command = "help"
      # puts "command par défaut positionnée à #{@command}"

      # Soit on en affiche la liste et on quitte :
      # (après avoir affiché la valeur des options selon demande)
      show_opts() if @level >= 2 || @gopts.show_opts
      puts "Il manque la commande à exécuter :".red
      # Optimist.educate # Ceci affiche la syntaxe complète du programme
      Args.show_commands
      puts "\n(TODO diam : remplacer \"optimist\" par \"optparse\" dans le code)\n"

      exit(1)
    end

    if ! @@commands.keys.include?(@command)
      puts "Erreur : commande \"#{@command}\" non autorisée :".red
      Args.show_commands
      exit(1)
    end

    # On positionne l'attribut @funtestsfile
    # Recherche et vérification du fichier de test
    @funtestsfile = @gopts.funtestsfile
    if @funtestsfile == "UNDEF"
      # Fichier funtestfile non précisé : on prend son chemin par défaut
      @funtestsfile = "Funtests.json"
      if !File.exist?(@funtestsfile) && @command != "proto"
        # Pas trouvé : on cheche dans le sous-répetoire test
        if File.exist?("test/Funtests.json")
          @funtestsfile = "test/Funtests.json"
        end
      end
    end

    # Si @funtestsfile est un répertoire, on y ajoute le fichier par défaut
    if File.directory?(@funtestsfile)
        @funtestsfile += "/Funtests.json"
    end

    # ici, @funtestsfile ne peut pas être UNDEF
    # le fichier @funtestsfile ne doit exister que pour certaines commmandes.
    # if  !File.exist?(@funtestsfile) && @command != "proto"
    if  !File.exist?(@funtestsfile) && %w(list run).include?(@command)
        msg="fichier de test inexistant ou illisible : \"#{@funtestsfile}\"."
        Optimist.die :funtestsfile, msg
    end

    # Exemple de vérification d'un type String énuméré
    # levels = %w(debug info warn error fatal)
    # if ! levels.include?(@gopts.loglevel)
    #   Optimist.die :loglevel, "doit-être dans #{levels}"
    # end

    # On extrait les options spécifiques à la sous-commande
    @copts = self.send("command_#{@command}")

    # Vérif et exploitation du level (verbosité utilisateur)
    # Attention @level ne doit-être écrasé que si copts.level* existe
    if @copts.level3
      @level = 3
    elsif @copts.level2
      @level = 2
    elsif @copts.level
      @level = @copts.level
    end
    @testpats = ARGV
  end

  def command_help
    # Traitement immédiat (pas d'options à gérer)
    puts usage
    exit
  end

  def command_list
    cmd = @command # pour passer dans la closure
    opts = Optimist::options do
      # l'attribut @command ne serait pas accessible ici (nil) !
      banner "options for #{cmd}"
      opt :level, "Niveau de verbosité" , :short => "-L", :default => 1
      opt :level2, "Niveau de verbosité à 2" , :short => "-v"
      opt :level3, "Niveau de verbosité à 3" , :short => "-V"
    end
  end

  def command_run
    cmd = @command # pour passer dans la closure
    opts = Optimist::options do
      # l'attribut @command ne serait pas accessible ici (nil) !
      banner "options for #{cmd}"
      opt :dry_run, "no execution", :short => "-n", :type => :boolean, :default => false
      opt :level, "Niveau de verbosité" , :short => "-L", :default => 1
      opt :level2, "Niveau de verbosité à 2" , :short => "-v"
      opt :level3, "Niveau de verbosité à 3" , :short => "-V"
    end
  end

  def command_proto
    cmd = @command # pour passer dans la closure
    opts = Optimist::options do
      # l'attribut @command ne serait pas accessible ici (nil) !
      banner "options for #{cmd}"
      opt :force, "ecrase existing file (- pour stdout)",
          :short => "-f", :type => :boolean, :default => false
      opt :format, "format d'affichage (si stdout) json ou hash",
          :short => "-F", :default => "json"
    end

    # S'il y a un argument supplémentaire : on prend le premier comme le nom du
    # fichier @funtestsfile à crééer
    if ARGV.length >= 1
      @funtestsfile = ARGV.shift
    end

    if !%w(json hash).include?(opts.format)
      Optimist.die :format, "#{opts.format} : seules les valeurs json et hash son autorisées"
    end

    if opts.format != "json" && @funtestsfile != "-"
      Optimist.die :format, "#{opts.format} : seules la valeur json est autorisée pour créer un fichier json"
    end
    return opts
  end

  # retourne le nom de l'action complèe à partir de son abréviation
  # Si ce n'est pas une abréviation possible, retourne la chaine initiale
  #
  def unabbrev(abbrev)
    actions = @@commands.keys
    abbrevs = actions.abbrev # hash de toutes les abbrevs possibles
    if abbrevs.has_key?(abbrev)
      return abbrevs[abbrev]
    else
      # return nil
      return abbrev
    end
  end

  def show_opts
    puts "Valeur de options globales :"
    puts @gopts.ai
    puts "Valeur de options pour la commande actives #{@command.ai} :"
    puts @copts.ai
    puts "valeur des attributs de Args : #{@level}"
    puts "  testpats = #{@testpats.ai}"
    puts "  funtestsfile = #{@funtestsfile.ai}"
    puts "  level = #{@level}"
  end
end

def usage
  # Attention notation "~" récente : Ruby 2.2+ (?)
  txt = <<~EOT
    ft=FunTest: Test fonctionnel d'une application à partir d'un fichier Funtests.json
    SYNTAXE : #{$EXE_NAME} -h # pour en savoir plus
    DESCRIPTION :
        Cette commande facilite les tests fonctionnels réalisés par l'appel
        d'une commande unix.
        Elle permet de lister ou d'exécuter une selection de tests.
        Chaque test est décrit dans un fichier (par défaut Funtest.json situé
        dans le répertoire courant).

    EXEMPLE :
        #{$EXE_NAME} proto
        => recopie un modèle de Funtests.json dans le répertoire courant
        #{$EXE_NAME} proto Funtests_dist.json
        => idem mais sous le nom de Funtests_dist.json
        #{$EXE_NAME} proto -
        => affichage du proto json sur la sortie standard
        #{$EXE_NAME} proto - --hash
        => idem mais après conversion en dictionnaire (ty Hash)
        #{$EXE_NAME} ou #{$EXE_NAME} -h
        => affiche la syntaxe
        #{$EXE_NAME} list
        => liste des tests disponibles dans le fichier json
        #{$EXE_NAME} list -v (-V)
        => liste détaillée (très détaillée) des tests disponibles dans le fichier json
        #{$EXE_NAME} list ls
        => liste détaillée des tests contenant la chaine "ls" dans le sid ou la desc
        #{$EXE_NAME} -f mytests/Funtests_bis.json list
        => utiliser le fichier Funtests_bis.json du sous-répertoire mytests
        #{$EXE_NAME} -f .. list
        => utiliser le fichier Funtests.json du répertoire parent
        #{$EXE_NAME} run  cplex -n     (ou --dry-run)
        => exécute "à blanc" les tests contenant cplex dans le sid ou dans desc
        #{$EXE_NAME} run  cplex -L 4   (ou --level)
        => exécute les tests et affiche leur exécution même si pas d'erreur

    CRÉATION D'UN FICHIER DE Funtests. POUR UNE APPLI UNIX QUELCONQUE
      - créer un fichier Funtests. avec la commande
        #{$EXE_NAME} proto
      - compléter le fichier Funtests.json avec vos propres règles.
      - le plus dur est de définir l'expression régulière qui valide le
        résultat de la commande testée.
          https://www.regular-expressions.info/quickstart.html
          https://www.regular-expressions.info/modifiers.html

    TODO :
       - s'il y a au moins un echec, enregistrer les détails dans un fichier log daté
       - créer attribut json exp_rpats au lieu de exp_rpat pour passer **une liste**
         de patternes à sélectionner
       - prévoir option --exclude (alias -x) pour ignorer des patternes de test
         (appelable de multiple fois), par exemple pour ignorer des tests trop longs.
       - compléter pour pouvoir tester des fichiers générés par une commande
         - prévoir de tester le nom d'un fichier généré
         - prévoir de tester le contenu d'un fichier généré
       - à exploiter pour valider une nouvelle installation de Julia, CPLEX, ...
  EOT
  return txt
end

# def myrun(help_string, result_string, cmd, dry_run=false )
def myrun(test, args )

  dry_run = args.copts.dry_run

  puts ""
  puts "=========================================================="
  if args.level >= 2
    puts test.to_s_long(dump: true)
    puts ""
  else
    puts "START: test #{test.to_s}".blue.underline
    puts test.cmd.bold
  end

  res_txt = ""
  if args.level >= 4
    loglevel = :info # pour afficher la sortie de la commande
  else
    loglevel = :warn
  end
  opts = {
    :loglevel => loglevel,
    :dry_run => dry_run,
    :onerror => :logerror,
    :prefix => "RUN: ",
    :redpat => /error|erreur|warning/i,  # coloriera les erreur/warning en rouge
    :capturevar => res_txt,
  }
  #
  # LA COMMANDE QUI FAIT LE BOULOT !
  #
  starttime = Time.now
  Madi.run_command(test.cmd, opts)
  res_duration = (Time.now - starttime).round(2)

  puts "=========================================================="

  if dry_run
    return
  end

  ###################################################################
  ###################################################################
  ###################################################################
  # Construction de la chaine d'affichage finale (pour tous les test
  # d'une commande)
  msgs = ""

  rpat_ok = true
  length_ok = true
  duration_ok = true

  ###################################################################
  # Validation de la pattern
  if test.exp_rpat
    rpat_ok = res_txt.match?(test.exp_rpat)
    msg = "=== exp_rpat :   \"#{test.exp_rpat}\" => "
    if rpat_ok
        msg << "OK "
        msgs << msg.green << "\n"
    else
        msg << "ECHEC "
        msgs << msg.red.on_light_yellow.bold << "\n"
    end
  end

  ###################################################################
  # Validation de la longueur du résultat
  if test.exp_length != -1
    msg = "=== exp_length=#{test.exp_length} (margin=#{test.exp_length_margin}) " <<
          "result:  #{res_txt.length} => "
    if test.exp_length_margin != -1
      val_min = test.exp_length - test.exp_length_margin
      val_max = test.exp_length + test.exp_length_margin
      length_ok = res_txt.length.between?(val_min, val_max)
    else
    end

    if length_ok
        msg << "OK "
        msgs << msg.green << "\n"
    else
        msg << "ECHEC "
        msgs << msg.red.on_light_yellow.bold << "\n"
    end
  end

  ###################################################################
  # Validation de la durée d'exécution
  if test.exp_duration != -1
    msg = "=== exp_duration=#{test.exp_duration}s (margin=#{test.exp_duration_margin}s) " <<
          "result:  #{res_duration}s => "
    if test.exp_duration_margin != -1
      val_min = test.exp_duration - test.exp_duration_margin
      val_max = test.exp_duration + test.exp_duration_margin
      duration_ok = res_duration.between?(val_min, val_max)
    end

    if duration_ok
        msg << "OK "
        msgs << msg.green << "\n"
    else
        msg << "ECHEC "
        msgs << msg.red.on_light_yellow.bold << "\n"
    end
  end

  msgs <<  "==========================================================\n"

  #############################################################################
  # Rappel éventuel du test effectué (fait si la sortie du test est affichée)
  if !rpat_ok || !length_ok
      puts "==== BEGIN RÉSULTAT ERRONÉ : =====".red
      puts res_txt
      puts "==== END RÉSULTAT ERRONÉ. =====".red
      puts "RAPPEL DU TEST ÉCHOUÉ :".red
      puts test.to_s_long
  elsif args.level >= 4
      puts ""
      puts "RAPPEL: ----------------------------------------------------------"
      puts test.to_s_long
  end

  $results[test] = rpat_ok && length_ok && duration_ok
  puts msgs if args.level >= 1
end

###############################################################################
# BEGIN : Module indépendant de cette application
module Madi
  # run_command : extrait de run_mip_gotic.rb
  # exécute une ligne de commande unix et affiche le résultat ligne par ligne.
  #
  # Exemple
  #   Madi::run_command("ssh 147.250.33.10 ls -al",
  #       # :logger => nil,
  #       # :dry_run => true,
  #       :level => :warn,
  #       :onerror => :logerror,
  #       :prefix => "LL: ",
  #       :redpat => /^d/, # colorie les répertoires en rouge
  #   )
  # TODO
  # - prévoir une option callback pour traiter chaque ligne de sortie (e.g pour
  #   l'afficher ou détecter une erreur)
  #
  # HIST
  # - 12/10/2017 : création à partir de mon run_mip_gotic.rb (de 2011)
  # - 13/05/2020 : ajout capturevar (pour madi_clitest.rb)
  #
  def self.run_command(command, user_opts=Hash.new)
    def_opts = {
      :logger => Logger.new(STDOUT), # objet potentiellement inutilisé !
      :loglevel => :info,
      :dry_run => false,
      :onerror => :logerror,
      :prefix => "RUN: ",
      :redpat => /error|erreur|warning/i,  # coloriera les erreur/warning en rouge
      :capturevar => nil,
    }
    opts = OpenStruct.new(def_opts.merge user_opts)
    opts.logger.level = opts.loglevel

    # puts "#{opts.prefix}run_command : START:\n   #{command}" if opts.logger.info?
    puts "run_command : START:\n   #{command}" if opts.logger.info?
    if opts.dry_run
      puts "\n*** mode dry_run (commande suivante non exécutée).***"
      puts "*** #{command} ***"
      return true
    end

    # Manuel ruby
    # merge standard output and standard error using
    # spawn option.  See the document of Kernel.spawn.
    # IO.popen(["ls", "/", :err=>[:child, :out]]) {|ls_io|
    #     ls_result_with_error = ls_io.read
    # }
    #
    # On n'est pas en mode dry_run => exécution réelle
    IO.popen("#{command} 2>&1", "r+") do |io|
      while (line = io.gets)
        # ap line
        if opts.logger.info?
          linemsg = opts.prefix + line.chomp
          if line =~ opts.redpat
            # linemsg = "\033[32m#{linemsg}\033[0m"
            linemsg = "\033[31m#{linemsg}\033[0m"
            # linemsg = "\033[31m\033[43m#{linemsg}\033[0m"
          end
          puts linemsg
        end
        if opts.capturevar
          opts.capturevar << line
        else
          puts "opts.capturevar est nil !"
        end
      end
    end
    # $? est positionné par la sortie de IO.open
    puts "Retour de IO.popen : #{$?.ai}" if opts.logger.info?
    ok = $?.success?
    if ok
      puts "run_command : succes #{command}" if opts.logger.info?
    else
      case opts[:onerror].to_sym
      when :logerror
        # opts.logger.error "run_command : echec de #{command}"
        errmsg = "run_command : echec de #{command}"
        puts "\033[31m#{errmsg}\033[0m" if opts.logger.error
      when :raise
        raise "Erreur lors le l'exécution de #{command}"
      else
        raise "Parametre opts.onerror='#{opts.onerror}' non reconnu"
      end
    end
    return ok
  end
end # module Madi
# END : Module indépendant de cette application
###############################################################################

def main_list args
  tests = get_funtests(args.funtestsfile, testpats: args.testpats)
  if tests.empty?
    puts "Aucune tests à effectuer".red  if args.level >= 1
    return
  end

  # On calcul la longueur max du sid pour améliorer formatage
  sid_len = tests.max{|t1,t2| t1.sid.length <=> t2.sid.length}.sid.length
  tests.each_index do |id|
    test = tests[id]
    if args.level >= 2
      puts test.to_s_long(dump: (args.level >= 3))
      # dump = args.level >= 3
      # puts test.to_s_long(dump: dump)
    else
    #   puts  "@%- #{sid_len}s: %s" % [test.sid, test.desc]
      print  "@%- #{sid_len}s: ".blue % [test.sid]
      puts  test.desc
    end
  end
end

def main_run args
  tests = get_funtests(args.funtestsfile, testpats: args.testpats)
  if tests.empty?
    puts "Aucune tests à effectuer".red if args.level >= 1
    return
  end

  # TODO : mémoriser le résultat pour affiher une synthèse finale.
  $results = Hash.new
  for test in tests
    $results[test] = false # passera éventuellemenet à triu en fin de test
    myrun test, args
  end
  # Affichage de la synthèse des test
  puts "Synthèse des tests"
  nb_oks = 0
  for test in tests
    ok = $results[test]
    nb_oks += 1 if ok
    if ok
      puts "ok: #{test.to_s}".green
    else
      puts "ko: #{test.to_s}".red
    end
  end
  nb_kos = tests.length - nb_oks
  print "Nombre de test effectués #{tests.length} "
  print " réussis: #{nb_oks} ".green
  if nb_kos == 0
    puts " echecs: #{nb_kos} ".green
  else
    puts " echecs: #{nb_kos} ".red
  end
end

def main_proto args
  # Recherche du fichier "Funtests.json" de la distribution
  appdir = File.dirname(File.realpath(__FILE__))
  src_path = appdir + "/Funtests.json"
  if !File.exist?(src_path)
    puts "Erreur : le fichier Funtests.json prototype est introuvable.".red
    puts "Il devrait être dans le même répertoire que #{$EXE_NAME}.".red
    puts "Est-il bien installé ?".red
    exit 1
  end

  # vérif de l'écraseent éventuel du fichier de destination
  dst_path = args.funtestsfile

  if dst_path == "-"
    json_txt = File.read(src_path)
    if args.copts.format == "json"
      puts json_txt
    elsif args.copts.format == "hash"
      specs = JSON.parse(json_txt)
      puts specs.ai
    end
  elsif !File.exist?(dst_path)
    FileUtils.cp src_path, dst_path, preserve:true, verbose:true
  elsif args.copts.force
    puts "Écrasement du fichier \"#{dst_path}\" !" if args.level >= 1
    FileUtils.cp src_path, dst_path, preserve:true, verbose:true
  else
    puts "Le fichier existe déjà \"#{dst_path}\" !".red
    puts "Ajouter l'option --force pour l'écraser".red
  end
end

def main args
  date1 = Time.now
  args = Args.new()
  args.show_opts if args.gopts.show_opts

  for tpat in args.testpats
    puts "test_pat : #{tpat}"
  end

  # On sous-traite au main de la commande (e.g. main_list)
  self.send("main_#{args.command}", args)
  date2 = Time.now
  puts "Durée totale des tests : " + (date2-date1).round(2).to_s + "s"

  exit
end

main ARGV
