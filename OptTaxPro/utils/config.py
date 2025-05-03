import configparser

def add_arg(config, argv, section, keys, v_type="str"):
	def __add(config, argv, section, key, v_type):
		if not config.has_option(section, key):
			return False

		value = (
#			config[section].getint(key) if v_type == 'int' else \
#			config[section].getfloat(key) if v_type == 'float' else \
			config[section].getboolean(key) if v_type == 'bool' else
			config[section].get(key)
		)
		if v_type == 'list':
			value = eval(value)
			argv += [key] + [str(v) for v in value]
		elif v_type == 'bool':
			if value:
				argv += [key]
		else:
			argv += [key, value]

		return argv

	if isinstance(keys, list):
		for key in keys:
			ret = __add(config, argv, section, key, v_type)
			if isinstance(ret, list):
				break
	else:
		ret = __add(config, argv, section, keys, v_type)

	return ret if isinstance(ret, list) else argv

def _parse_preprocess(config, argv):
	argv = add_arg(config, argv, 'PREPROCESS', ['-m', '--min-len'], v_type='int')
	argv = add_arg(config, argv, 'PREPROCESS', ['-M', '--max-len'], v_type='int')
	argv = add_arg(config, argv, 'PREPROCESS', '-O', v_type='int')
	argv = add_arg(config, argv, 'PREPROCESS', '-e', v_type='float')

	argv = add_arg(config, argv, 'PREPROCESS', '--qc-cutoff', v_type='float')

	argv = add_arg(config, argv, 'PREPROCESS', '--single-cutoff', v_type='float')
	argv = add_arg(config, argv, 'PREPROCESS', '--single-threads', v_type='int')

	return argv

def _parse_cluster(config, argv):
	argv = add_arg(config, argv, 'CLUSTER', '--otu-cutoff', v_type='float')
	argv = add_arg(config, argv, 'CLUSTER', '--otu-threads', v_type='int')

	return argv

def _parse_classify(config, argv):
	argv = add_arg(config, argv, 'CLASSIFY', ['-d', '--db'], v_type='str')
	argv = add_arg(config, argv, 'CLASSIFY', '--search-threads', v_type='int')
	argv = add_arg(config, argv, 'CLASSIFY', '--search-cutoffs', v_type='list')

	argv = add_arg(config, argv, 'CLASSIFY', ['-r', '--ranks'], v_type='list')
	argv = add_arg(config, argv, 'CLASSIFY', '--assign-cutoffs', v_type='list')
	argv = add_arg(config, argv, 'CLASSIFY', '--acc2taxid', v_type='str')
	argv = add_arg(config, argv, 'CLASSIFY', '--hsg', v_type='str')
	argv = add_arg(config, argv, 'CLASSIFY', '--remove-self', v_type='bool')
	argv = add_arg(config, argv, 'CLASSIFY', '--add-name', v_type='bool')

	argv = add_arg(config, argv, 'CLASSIFY', '--u_dir', v_type='str')
	argv = add_arg(config, argv, 'CLASSIFY', '--uc-suffix', v_type='str')

	return argv

def _parse_profile(config, argv):
	argv = add_arg(config, argv, 'PROFILE', '--profile-ranks', v_type='list')
	argv = add_arg(config, argv, 'PROFILE', '--base-col', v_type='str')
	argv = add_arg(config, argv, 'PROFILE', '--filtering-cutoffs', v_type='list')
	argv = add_arg(config, argv, 'PROFILE', '--filtering-pivot', v_type='str')
	argv = add_arg(config, argv, 'PROFILE', '--output-prefix', v_type='str')

	argv = add_arg(config, argv, 'PROFILE', '--u_dir', v_type='str')
	argv = add_arg(config, argv, 'PROFILE', '--uc-suffix', v_type='str')

	return argv

def parse_cfg(file_name, prog):
	config = configparser.ConfigParser()
	config.read(file_name)

	argv = [prog]
	argv = add_arg(config, argv, 'DEFAULT', ['-h', '--help'], v_type='bool')
	argv = add_arg(config, argv, 'DEFAULT', ['-q', '--quiet'], v_type='bool')
	argv = add_arg(config, argv, 'DEFAULT', '--log', v_type='str')
	argv = add_arg(config, argv, 'DEFAULT', ['-p', '--processes'], v_type='int')
	argv = add_arg(config, argv, 'DEFAULT', '--suffix', v_type='str')
	argv = add_arg(config, argv, 'DEFAULT', '--t_dir', v_type='str')
	argv = add_arg(config, argv, 'DEFAULT', '--o_dir', v_type='str')
	argv = add_arg(config, argv, 'DEFAULT', '--i_dir', v_type='str')

	if prog == 'preprocess':
		argv = _parse_preprocess(config, argv)
	elif prog == 'cluster':
		argv = _parse_cluster(config, argv)
	elif prog == 'classify':
		argv = _parse_classify(config, argv)
	elif prog == 'profile':
		argv = _parse_profile(config, argv)
	elif prog == 'alltheway':
		argv = _parse_preprocess(config, argv)
		argv = _parse_cluster(config, argv)
		argv = _parse_classify(config, argv)
		argv = _parse_profile(config, argv)

	return argv
