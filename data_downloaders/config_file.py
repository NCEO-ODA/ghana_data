#!/usr/bin/env python
"""Add a configuration file to click command line parser
See here: https://serge-m.github.io/click-config-parsers.html
"""
import click
import yaml



class CustomCommandClass(click.Command):

    def invoke(self, ctx):
        config_file = ctx.params[config_file_param_name]
        if config_file is not None:
            with open(config_file) as f:
                config_data = yaml.load(f)
                for param, value in ctx.params.items():
                    if value is None and param in config_data:
                        ctx.params[param] = config_data[param]

        return super(CustomCommandClass, self).invoke(ctx)

return CustomCommandClass