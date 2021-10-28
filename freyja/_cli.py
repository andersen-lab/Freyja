import click
#from freyja.convert_paths2barcodes import conv_to_barcodes


@click.group()
def cli():
    pass


@cli.command()
@click.option('--input-csv', required=True, type=click.Path(exists=True),
              help="This is the magical input data")
@click.option('--output-csv', required=True, type=click.Path(exists=False),
              help="This is the magical output")
def convert_paths_to_barcodes(input_csv, output_csv):
    # load our data frame#
    df = None
    result = conv_to_barcodes(df)
    result.to_csv(output_csv)


if __name__ == '__main__':
    cli()
