class DecayChannelNotFoundError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Could not determine the decay channel') 


def get_decay_channel_tag(args):
    if args.combination:
        tag = 'combination'
    elif args.hgg:
        tag = 'hgg'
    elif args.hzz:
        tag = 'hzz'
    elif args.hbb:
        tag = 'hbb'
    elif args.combWithHbb:
        tag = 'combWithHbb'
    else:
        raise DecayChannelNotFoundError()
    return tag
